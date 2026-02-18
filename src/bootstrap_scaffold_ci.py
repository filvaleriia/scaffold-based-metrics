"""
Bootstrap confidence intervals (CI) for scaffold-based metrics.

This module is a cleaned-up and configurable version of the original script you provided.
Key features:
- Two independent parallelism knobs:
  1) job_workers: number of parallel (receptor, split, scaffold, generator) jobs
  2) scaffold_workers: number of processes used for SMILES -> scaffold conversion inside each job
- Per-input-file scaffold caching stored as NumPy .npy (dtype=object)
- Per-job outputs stored under: results_dir/<receptor>/<scaffold>_scaffolds/<split>/<generator>/
- Optional saving of bootstrap distributions (needed for significance testing)

Metrics computed (per bootstrap replicate, after averaging across available clusters):
- SED: unique scaffolds / n_subsample
- ASER: recall matches (with replacement) / n_subsample
- RS: unique active scaffolds / unique recall scaffolds

Expected file layout under data_folder (same as your original):
- Output sets:
  data/output_sets/<receptor>/<generator>/cOS_<generator>_<split>_<cluster>_one_column.csv
- Recall sets:
  data/input_recall_sets/<receptor>/cRS_<receptor>_<split>_<cluster>.csv
"""

from __future__ import annotations

import os
import time
import argparse
from dataclasses import dataclass
from pathlib import Path
from functools import partial
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles, MakeScaffoldGeneric


ScaffoldType = Literal["csk", "murcko"]


# =========================
# Configuration dataclass
# =========================

@dataclass(frozen=True)
class BootstrapConfig:
    clusters: list[int]
    n_subsample: int
    n_bootstrap: int
    alpha: float

    results_dir: str
    data_folder: str

    job_workers: int
    scaffold_workers: int

    chunksize: int = 2000
    use_cache: bool = True
    show_progress: bool = True

    save_bootstrap_samples: bool = True  # saves per-metric bootstrap arrays per job


# =========================
# SMILES -> scaffold
# =========================

def convert_to_scaffold(smiles: str, scaffold_type: ScaffoldType) -> str | None:
    """Convert a SMILES to a scaffold SMILES. Returns None if conversion fails."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if scaffold_type == "csk":
            sc = MurckoScaffoldSmiles(Chem.MolToSmiles(MakeScaffoldGeneric(mol)))
        elif scaffold_type == "murcko":
            sc = MurckoScaffoldSmiles(Chem.MolToSmiles(mol))
        else:
            return None
        return sc if sc else None
    except Exception:
        return None


def _cache_path(smiles_file: str, scaffold_type: ScaffoldType) -> str:
    """Cache path: <file>.<scaffold>.scaffolds.npy"""
    p = Path(smiles_file)
    return str(p.with_suffix(p.suffix + f".{scaffold_type}.scaffolds.npy"))


def _read_smiles(filepath: str) -> list[str]:
    with open(filepath) as f:
        return f.read().splitlines()


def smiles_to_scaffolds(
    filepath: str,
    scaffold_type: ScaffoldType,
    ncpus: int,
    chunksize: int,
    use_cache: bool,
    show_progress: bool,
    label: str = "",
) -> np.ndarray:
    """
    Convert a SMILES file (one SMILES per line) to an array of scaffolds.
    Uses multiprocessing Pool for parallel conversion and stores a .npy cache.
    """
    cache = _cache_path(filepath, scaffold_type)

    if use_cache and os.path.exists(cache):
        arr = np.load(cache, allow_pickle=True)
        return arr

    smiles = _read_smiles(filepath)
    n = len(smiles)
    if n == 0:
        return np.asarray([], dtype=object)

    worker = partial(convert_to_scaffold, scaffold_type=scaffold_type)

    scaffolds: list[str] = []
    with Pool(processes=ncpus) as pool:
        it = pool.imap_unordered(worker, smiles, chunksize=chunksize)
        if show_progress:
            it = tqdm(it, total=n, unit="mol", desc=f"{label} scaffolds[{scaffold_type}]", smoothing=0.05)
        for sc in it:
            if sc is not None:
                scaffolds.append(sc)

    arr = np.asarray(scaffolds, dtype=object)
    if use_cache:
        np.save(cache, arr, allow_pickle=True)
    return arr


# =========================
# Bootstrap per cluster
# =========================

def bootstrap_cluster(
    output_scaffolds: np.ndarray,
    recall_unique: set[str],
    n_recall_unique: int,
    n_sub: int,
    B: int,
    seed: int,
) -> dict[str, np.ndarray]:
    """
    Bootstrap with replacement from output_scaffolds.
    Returns arrays of length B for metrics: RS, SED, ASER.
    """
    rng = np.random.default_rng(seed)

    if len(output_scaffolds) == 0:
        z = np.zeros(B, dtype=float)
        return {"RS": z.copy(), "SED": z.copy(), "ASER": z.copy()}

    idx = rng.integers(0, len(output_scaffolds), size=(B, n_sub))
    samp = output_scaffolds[idx]  # (B, n_sub), dtype=object

    T = np.zeros(B, dtype=float)
    S = np.zeros(B, dtype=float)
    A = np.zeros(B, dtype=float)

    # NOTE: this loop is intentionally explicit for clarity and easy debugging.
    for b in range(B):
        s = samp[b]
        uniq = set(s.tolist())
        n_unique_out = len(uniq)

        n_matches = sum(x in recall_unique for x in s)     # with replacement
        n_unique_active_out = sum(x in recall_unique for x in uniq)

        S[b] = n_unique_out / n_sub
        A[b] = n_matches / n_sub
        T[b] = (n_unique_active_out / n_recall_unique) if n_recall_unique > 0 else 0.0

    return {"RS": T, "SED": S, "ASER": A}


# =========================
# IO helpers
# =========================

def job_output_dir(results_dir: str, receptor: str, split: str, scaffold: str, generator: str) -> Path:
    return Path(results_dir) / receptor / f"{scaffold}_scaffolds" / split / generator


def save_job_outputs(
    results_dir: str,
    receptor: str,
    split: str,
    scaffold: str,
    generator: str,
    used_clusters: list[int],
    df_job: pd.DataFrame,
    bootstrap_samples: dict[str, np.ndarray] | None,
    cfg: BootstrapConfig,
) -> None:
    out_dir = job_output_dir(results_dir, receptor, split, scaffold, generator)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_job.to_csv(out_dir / "bootstrap_ci.csv", index=False)

    # Optional: store bootstrap arrays for significance testing (one per metric)
    if bootstrap_samples is not None:
        bs_dir = out_dir / "bootstrap_samples"
        bs_dir.mkdir(exist_ok=True)
        for metric, arr in bootstrap_samples.items():
            np.save(bs_dir / f"{metric}.npy", np.asarray(arr, dtype=float))

    # Human-readable metadata
    meta = out_dir / "meta.txt"
    with open(meta, "w") as f:
        f.write(f"receptor={receptor}\n")
        f.write(f"split={split}\n")
        f.write(f"scaffold={scaffold}\n")
        f.write(f"generator={generator}\n")
        f.write(f"clusters_used={used_clusters}\n")
        f.write(f"n_bootstrap={cfg.n_bootstrap}\n")
        f.write(f"n_subsample={cfg.n_subsample}\n")
        f.write(f"alpha={cfg.alpha}\n")
        f.write(f"job_workers={cfg.job_workers}\n")
        f.write(f"scaffold_workers={cfg.scaffold_workers}\n")
        f.write(f"chunksize={cfg.chunksize}\n")
        f.write(f"use_cache={cfg.use_cache}\n")
        f.write(f"save_bootstrap_samples={cfg.save_bootstrap_samples}\n")


# =========================
# One job
# =========================

def _paths_for_cluster(data_folder: str, receptor: str, split: str, generator: str, cluster: int) -> tuple[str, str]:
    out_f = f"{data_folder}/data/output_sets/{receptor}/{generator}/cOS_{generator}_{split}_{cluster}_one_column.csv"
    rec_f = f"{data_folder}/data/input_recall_sets/{receptor}/cRS_{receptor}_{split}_{cluster}.csv"
    return out_f, rec_f


def run_job(
    receptor: str,
    split: str,
    scaffold: ScaffoldType,
    generator: str,
    cfg: BootstrapConfig,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    per_cluster: dict[str, list[np.ndarray]] = {m: [] for m in ["RS", "SED", "ASER"]}
    used_clusters: list[int] = []

    for c in cfg.clusters:
        out_f, rec_f = _paths_for_cluster(cfg.data_folder, receptor, split, generator, c)
        if not (os.path.exists(out_f) and os.path.exists(rec_f)):
            continue

        out_sc = smiles_to_scaffolds(
            out_f,
            scaffold_type=scaffold,
            ncpus=cfg.scaffold_workers,
            chunksize=cfg.chunksize,
            use_cache=cfg.use_cache,
            show_progress=cfg.show_progress,
            label=f"OS {generator} {split} c{c}",
        )
        rec_sc = smiles_to_scaffolds(
            rec_f,
            scaffold_type=scaffold,
            ncpus=cfg.scaffold_workers,
            chunksize=cfg.chunksize,
            use_cache=cfg.use_cache,
            show_progress=cfg.show_progress,
            label=f"RS {receptor} {split} c{c}",
        )

        rec_u = set(np.unique(rec_sc).tolist())
        n_recall_unique = len(rec_u)

        d = bootstrap_cluster(
            output_scaffolds=out_sc,
            recall_unique=rec_u,
            n_recall_unique=n_recall_unique,
            n_sub=cfg.n_subsample,
            B=cfg.n_bootstrap,
            seed=int(rng.integers(1_000_000_000)),
        )

        for m in per_cluster:
            per_cluster[m].append(d[m])
        used_clusters.append(c)

    # If no data, write an error row (keeps downstream processing simple)
    if not used_clusters:
        df_err = pd.DataFrame([{
            "receptor": receptor,
            "split": split,
            "scaffold": scaffold,
            "generator": generator,
            "metric": None,
            "mean": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "clusters": "",
            "error": "no clusters found (missing files)",
        }])
        save_job_outputs(
            cfg.results_dir, receptor, split, scaffold, generator,
            used_clusters=[], df_job=df_err, bootstrap_samples=None, cfg=cfg
        )
        return df_err

    # Aggregate across clusters: for each bootstrap replicate, average metric across clusters
    rows = []
    bootstrap_samples: dict[str, np.ndarray] | None = {} if cfg.save_bootstrap_samples else None

    for m in per_cluster:
        mat = np.vstack(per_cluster[m])      # (k_clusters, B)
        avg = mat.mean(axis=0)               # (B,)

        mean_v = float(avg.mean())
        lo = float(np.quantile(avg, cfg.alpha / 2))
        hi = float(np.quantile(avg, 1.0 - cfg.alpha / 2))

        rows.append({
            "receptor": receptor,
            "split": split,
            "scaffold": scaffold,
            "generator": generator,
            "metric": m,
            "mean": mean_v,
            "ci_low": lo,
            "ci_high": hi,
            "clusters": ",".join(map(str, used_clusters)),
            "error": None,
        })

        if bootstrap_samples is not None:
            bootstrap_samples[m] = avg.astype(float, copy=False)

    df_job = pd.DataFrame(rows)

    save_job_outputs(
        cfg.results_dir, receptor, split, scaffold, generator,
        used_clusters=used_clusters, df_job=df_job,
        bootstrap_samples=bootstrap_samples, cfg=cfg
    )
    return df_job


# =========================
# Orchestration
# =========================

def run_all_jobs(
    receptors: Iterable[str],
    splits: Iterable[str],
    scaffolds: Iterable[ScaffoldType],
    generators: Iterable[str],
    cfg: BootstrapConfig,
    base_seed: int = 42,
) -> pd.DataFrame:
    jobs = [(r, s, sc, g) for r in receptors for s in splits for sc in scaffolds for g in generators]

    t0 = time.time()
    out: list[pd.DataFrame] = []

    with ProcessPoolExecutor(max_workers=cfg.job_workers) as ex:
        futs = []
        for i, (r, s, sc, g) in enumerate(jobs):
            futs.append(ex.submit(run_job, r, s, sc, g, cfg, base_seed + i * 1000))

        for f in tqdm(as_completed(futs), total=len(futs), desc="Jobs", unit="job"):
            out.append(f.result())

    df = pd.concat(out, ignore_index=True)

    Path(cfg.results_dir).mkdir(parents=True, exist_ok=True)
    summary_path = Path(cfg.results_dir) / "bootstrap_ci_scaffold_metrics_summary.csv"
    df.to_csv(summary_path, index=False)

    dt = time.time() - t0
    print(f"[DONE] Wrote: {summary_path} ({dt/60:.2f} min)")
    return df


# =========================
# CLI
# =========================

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Bootstrap CI for scaffold metrics (configurable parallelism).")

    p.add_argument("--data-folder", default="../../", help="Base folder containing the 'data/' directory.")
    p.add_argument("--results-dir", default="results_bootstrap_ci", help="Where to save outputs.")

    p.add_argument("--job-workers", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                   help="Parallel jobs in ProcessPoolExecutor.")
    p.add_argument("--scaffold-workers", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                   help="Processes for SMILES->scaffold conversion inside a job.")

    p.add_argument("--n-subsample", type=int, default=250_000, help="Subsample size per bootstrap replicate.")
    p.add_argument("--n-bootstrap", type=int, default=300, help="Number of bootstrap replicates.")
    p.add_argument("--alpha", type=float, default=0.05, help="CI level: (alpha/2, 1-alpha/2).")

    p.add_argument("--clusters", default="0,1,2,3,4", help="Comma-separated cluster IDs, e.g. 0,1,2,3,4.")
    p.add_argument("--chunksize", type=int, default=2000, help="multiprocessing imap_unordered chunksize.")
    p.add_argument("--no-cache", action="store_true", help="Disable scaffold caching.")
    p.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bars.")
    p.add_argument("--no-bootstrap-samples", action="store_true", help="Do not save bootstrap distributions per job.")
    return p


def main_cli() -> None:
    p = build_argparser()
    args = p.parse_args()

    clusters = [int(x.strip()) for x in args.clusters.split(",") if x.strip()]

    cfg = BootstrapConfig(
        clusters=clusters,
        n_subsample=int(args.n_subsample),
        n_bootstrap=int(args.n_bootstrap),
        alpha=float(args.alpha),
        results_dir=str(args.results_dir),
        data_folder=str(args.data_folder),
        job_workers=int(args.job_workers),
        scaffold_workers=int(args.scaffold_workers),
        chunksize=int(args.chunksize),
        use_cache=not bool(args.no_cache),
        show_progress=not bool(args.no_progress),
        save_bootstrap_samples=not bool(args.no_bootstrap_samples),
    )

    receptors = ["Glucocorticoid_receptor", "Leukocyte_elastase"]
    splits = ["dis", "sim"]
    scaffolds: list[ScaffoldType] = ["csk", "murcko"]
    generators = [
        "Molpher",
        "DrugEx_GT_epsilon_0.1",
        "DrugEx_GT_epsilon_0.6",
        "DrugEx_RNN_epsilon_0.1",
        "DrugEx_RNN_epsilon_0.6",
        "GB_GA_log_p_mut_r_0.01",
        "GB_GA_log_p_mut_r_0.5",
        "GB_GA_mut_r_0.01",
        "GB_GA_mut_r_0.5",
        "REINVENT",
        "addcarbon",
    ]

    run_all_jobs(receptors, splits, scaffolds, generators, cfg)


if __name__ == "__main__":
    main_cli()
