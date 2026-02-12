"""
Molecular Scaffold Metrics Calculator

Calculates similarity metrics (RS, SED, ASER) between output and recall sets 
using Murcko scaffolds or CSK (core scaffold with generic atoms).

Usage:
    Command line:
        python metrics_scaffold.py --type_cluster dis --type_scaffold csk --generator addcarbon_250k --receptor Glucocorticoid_receptor
    
    Jupyter:
        from metrics_scaffold import Metrics_scaffold
        mt = Metrics_scaffold('dis', 'csk', 'addcarbon_250k', 'Glucocorticoid_receptor')
        mt.calculate()
"""

import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
from pathlib import Path
import argparse
from typing import List, Tuple, Optional

from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles, MakeScaffoldGeneric


def convert_to_scaffold(smiles_str: str, scaffold_type: str) -> Optional[str]:
    """
    Convert a SMILES string to its scaffold representation.
    
    Args:
        smiles_str: Input SMILES string
        scaffold_type: Type of scaffold ('csk' or 'murcko')
    
    Returns:
        Scaffold SMILES string or None if conversion fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return None
        
        if scaffold_type == 'csk':
            # Generic Murcko scaffold (core structure with generic atoms)
            scaffold = MurckoScaffoldSmiles(
                Chem.MolToSmiles(MakeScaffoldGeneric(mol))
            )
        elif scaffold_type == 'murcko':
            # Standard Murcko scaffold
            scaffold = MurckoScaffoldSmiles(Chem.MolToSmiles(mol))
        else:
            return None
        
        return scaffold if scaffold != '' else None
    except Exception:
        return None


def create_matching_dataframe(output_scaffolds: pd.DataFrame, recall_scaffolds: List[str]) -> pd.DataFrame:
    """
    Create matching dataframe showing how many times each recall scaffold appears in output.
    
    Args:
        output_scaffolds: DataFrame with output scaffolds
        recall_scaffolds: List of unique recall scaffolds
    
    Returns:
        DataFrame with columns: unique_scaffold_recall, count_of_occurrence, unique_indicator
    """
    result = pd.DataFrame({'scaffold': recall_scaffolds})
    counts = output_scaffolds[0].value_counts()
    
    result['count'] = result['scaffold'].apply(
        lambda x: counts[x] if x in counts.index else 0
    )
    result['unique_indicator'] = result['count'].apply(lambda x: 1 if x > 0 else 0)
    
    result.columns = ['unique_scaffold_recall', 'count_of_occurrence', 'unique_indicator']
    
    return result


class Metrics_scaffold:
    """
    Molecular scaffold metrics calculator.
    
    Attributes:
        type_cluster: Type of clustering ('dis' or 'sim')
        type_scaffold: Type of scaffold ('csk' or 'murcko')
        generator_name: Name of the molecule generator
        receptor: Receptor name
        data_folder: Base path to data folder (default: '')
        ncpus: Number of CPUs for multiprocessing (default: 1)
    """
    
    def __init__(
        self,
        type_cluster: str,
        type_scaffold: str,
        generator_name: str,
        receptor: str,
        data_folder: str = '',
        ncpus: int = 1
    ):
        self.type_cluster = type_cluster
        self.type_scaffold = type_scaffold
        self.generator_name = generator_name
        self.receptor = receptor
        self.data_folder = data_folder
        self.ncpus = ncpus
        
        # Will be populated during calculation
        self.number_of_calculation: Optional[int] = None
        self.output_set: Optional[pd.DataFrame] = None
        self.recall_set: Optional[pd.DataFrame] = None
        self.output_set_scaffolds: Optional[pd.DataFrame] = None
        self.recall_set_scaffolds: Optional[pd.DataFrame] = None
        self.unique_output_set: Optional[np.ndarray] = None
        self.unique_recall_set: Optional[np.ndarray] = None
        self.count_metrics: Optional[pd.DataFrame] = None
        self.results: pd.DataFrame = pd.DataFrame()

    def _get_output_dir(self) -> str:
        """Get output directory path."""
        return f"{self.data_folder}/data/results/{self.receptor}/{self.type_scaffold}_scaffolds/{self.type_cluster}/{self.generator_name}"

    def load(self, filepath_output_set: str, filepath_recall_set: str):
        """
        Load output and recall sets from files.
        
        Args:
            filepath_output_set: Path to output set file
            filepath_recall_set: Path to recall set file
        """
        #print(f"\nLoading data from:")
        #print(f"  Output: {filepath_output_set}")
        #print(f"  Recall: {filepath_recall_set}")
        
        with open(filepath_output_set, 'r') as f:
            output_set = f.read().splitlines()
        self.output_set = pd.DataFrame(output_set)
        
        with open(filepath_recall_set, 'r') as f:
            recall_set = f.read().splitlines()
        self.recall_set = pd.DataFrame(recall_set)
        
        #print(f"Loaded {len(self.output_set):,} output compounds")
        #print(f"Loaded {len(self.recall_set):,} recall compounds")

    def convert_to_scaffolds(self):
        """Convert SMILES to scaffolds using multiprocessing."""
        #print(f"\nConverting to {self.type_scaffold} scaffolds...")
        
        # Convert recall set
        #print("  Processing recall set...")
        with Pool(processes=self.ncpus) as pool:
            recall_results = pool.starmap(
                convert_to_scaffold,
                [(smiles, self.type_scaffold) for smiles in self.recall_set[0]]
            )
        self.recall_set_scaffolds = pd.DataFrame(recall_results).dropna()
        #print(f"    Recall scaffolds: {len(self.recall_set_scaffolds):,}")
        
        # Convert output set
        #print("  Processing output set...")
        with Pool(processes=self.ncpus) as pool:
            output_results = pool.starmap(
                convert_to_scaffold,
                [(smiles, self.type_scaffold) for smiles in self.output_set[0]]
            )
        self.output_set_scaffolds = pd.DataFrame(output_results).dropna()
        #print(f"    Output scaffolds: {len(self.output_set_scaffolds):,}")
        
        # Get unique scaffolds
        self.unique_output_set = self.output_set_scaffolds[0].unique()
        self.unique_recall_set = self.recall_set_scaffolds[0].unique()
        
        #print(f"  Unique output scaffolds: {len(self.unique_output_set):,}")
        #print(f"  Unique recall scaffolds: {len(self.unique_recall_set):,}")

    def calculate_metrics(self):
        """Calculate similarity metrics between output and recall sets."""
        #print("\nCalculating metrics...")
        
        # Convert to scaffolds
        self.convert_to_scaffolds()
        
        # Create matching dataframe
        self.count_metrics = create_matching_dataframe(
            self.output_set_scaffolds,
            self.unique_recall_set
        )
        
        # Calculate metrics
        USo = len(self.unique_output_set)  # Unique scaffolds in output
        SSo = len(self.output_set_scaffolds)  # Total scaffolds in output
        cASo = self.count_metrics['count_of_occurrence'].sum()  # Total compound matches
        
        # Unique active scaffolds in output (scaffolds that matched recall)
        try:
            UASo = self.count_metrics['unique_indicator'].value_counts()[1]
        except (KeyError, IndexError):
            UASo = 0
        
        UASr = len(self.count_metrics)  # Unique active scaffolds in recall
        
        # Calculate derived metrics
        try:
            RS = UASo / UASr
            RS_text = f"{UASo}/{UASr}"
        except ZeroDivisionError:
            RS = 0
            RS_text = f"0/{UASr}"
        
        SED = USo / SSo if SSo > 0 else 0
        ASER = cASo / SSo if SSo > 0 else 0
        
        # Create results dataframe
        self.results = pd.DataFrame({
            'name': [f"{self.generator_name}_{self.number_of_calculation}"],
            'type_cluster': [self.type_cluster],
            'scaffold': [self.type_scaffold],
            'SSo': [SSo],
            'RS_': [RS_text],
            'RS': [RS],
            'SED': [SED],
            'ASER': [ASER]
        })
        
        #print(f"\nMetrics calculated:")
        #print(f"  RS: {RS:.4f} ({RS_text})")
        #print(f"  SED:  {SED:.4f}")
        #print(f"  ASER:  {ASER:.4f}")
        
        return self.type_cluster, SSo, RS_text, RS, SED, ASER

    def save_results(self):
        """Save calculated metrics and scaffolds to files."""
        output_dir = self._get_output_dir()
        os.makedirs(output_dir, exist_ok=True)
        
        count_file = f"{output_dir}/count_of_occurrence_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv"
        metrics_file = f"{output_dir}/metrics_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv"
        output_scaffolds_file = f"{output_dir}/scaffolds_of_output_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv"
        recall_scaffolds_file = f"{output_dir}/scaffolds_of_recall_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv"
        
        self.count_metrics.to_csv(count_file, index=False)
        self.results.to_csv(metrics_file, index=False)
        self.output_set_scaffolds.to_csv(output_scaffolds_file, header=False, index=False)
        self.recall_set_scaffolds.to_csv(recall_scaffolds_file, header=False, index=False)
        
        #print(f"\nResults saved to: {output_dir}")

    def average_value(self, numbers: List[int]) -> pd.DataFrame:
        """
        Calculate average metrics across multiple clusters.
        
        Args:
            numbers: List of cluster numbers to average
        
        Returns:
            Combined DataFrame with mean values
        """
        output_dir = self._get_output_dir()
        
        # Load all cluster results
        file_paths = [
            f"{output_dir}/metrics_cluster_{num}_{self.type_cluster}_{self.generator_name}.csv"
            for num in numbers
        ]
        
        #print(f"\nAveraging results across {len(numbers)} clusters...")
        combined_df = pd.concat([pd.read_csv(path) for path in file_paths], ignore_index=True)
        
        # Calculate mean for numeric columns
        mean_values = combined_df.mean(numeric_only=True)
        
        # Create mean row
        mean_row = {
            'name': f"{self.generator_name}_mean",
            'type_cluster': self.type_cluster,
            'scaffold': self.type_scaffold,
            'SSo': mean_values.get('SSo', np.nan),
            'RS_': '-',
            'RS': mean_values.get('RS', np.nan),
            'SED': mean_values.get('SED', np.nan),
            'ASER': mean_values.get('ASER', np.nan)
        }
        
        # Append mean row
        combined_df = pd.concat([combined_df, pd.DataFrame([mean_row])], ignore_index=True)
        combined_df = combined_df.round(7)
        
        # Create formatted version for display
        formatted_df = combined_df.copy()
        for col in ['SSo']:
            formatted_df[col] = formatted_df[col].apply(
                lambda x: "{:,}".format(int(x)) if pd.notnull(x) else x
            )
        
        # Save results
        formatted_df.to_csv(f"{output_dir}/df_all_clusters_with_mean_with_coma.csv", index=False)
        combined_df.to_csv(f"{output_dir}/df_all_clusters_with_mean.csv", index=False)
        formatted_df.to_html(f"{output_dir}/df_all_clusters_with_mean.html", index=False)
        combined_df.tail(1).to_csv(f"{output_dir}/{self.generator_name}_mean_{self.type_scaffold}_{self.type_cluster}.csv", index=False)
        
        #print(f"Average results saved to: {output_dir}")
        
        return combined_df

    def calculate(self, cluster_range: range = range(5)):
        """
        Calculate metrics for all available clusters.
        
        Args:
            cluster_range: Range of cluster numbers to process (default: 0-4)
        
        Returns:
            DataFrame with results for all clusters including mean
        """
        print(f"\n{'='*60}")
        print(f"Starting calculation for {self.generator_name}")
        print(f"Receptor: {self.receptor}")
        print(f"Cluster type: {self.type_cluster}")
        print(f"Scaffold type: {self.type_scaffold}")
        print(f"CPUs: {self.ncpus}")
        print(f"{'='*60}")
        
        numbers = []
        
        for number in cluster_range:
            self.number_of_calculation = number
            
            output_file_path = (
                f"{self.data_folder}/data/output_sets/{self.receptor}/{self.generator_name}/"
                f"cOS_{self.generator_name}_{self.type_cluster}_{number}_one_column.csv"
            )
            recall_file_path = (
                f"{self.data_folder}/data/input_recall_sets/{self.receptor}/"
                f"cRS_{self.receptor}_{self.type_cluster}_{number}.csv"
            )
            
            if os.path.exists(output_file_path):
                print(f"\n--- Processing cluster {number} ---")
                self.load(output_file_path, recall_file_path)
                self.calculate_metrics()
                self.save_results()
                numbers.append(number)
            else:
                print(output_file_path)
                print(recall_file_path)
                print(f"\nSkipping cluster {number} (file not found)")
        
        if numbers:
            result = self.average_value(numbers)
            #print(f"\n{'='*60}")
            print(f"Calculation complete! Processed {len(numbers)} clusters")
            #print(f"{'='*60}")
            return result[['name', 'type_cluster', 'scaffold', 'RS', 'SED', 'ASER']].copy()
        else:
            #print('\nNo data found for calculation')
            return pd.DataFrame()


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description='Calculate molecular scaffold similarity metrics.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python metrics_scaffold.py --type_cluster dis --type_scaffold csk --generator addcarbon_250k --receptor Glucocorticoid_receptor
  python metrics_scaffold.py --type_cluster sim --type_scaffold murcko --generator Molpher_250k --receptor Glucocorticoid_receptor --ncpus 8
        """
    )
    
    parser.add_argument('--type_cluster', type=str, required=True,
                        help='Type of clustering (dis/sim)')
    parser.add_argument('--type_scaffold', type=str, required=True,
                        help='Type of scaffold (csk/murcko)')
    parser.add_argument('--generator', type=str, required=True,
                        help='Generator name (e.g., addcarbon_250k)')
    parser.add_argument('--receptor', type=str, required=True,
                        help='Receptor name (e.g., Glucocorticoid_receptor)')
    parser.add_argument('--data_folder', type=str, default='',
                        help='Base path to data folder (default: current directory)')
    parser.add_argument('--ncpus', type=int, default=1,
                        help='Number of CPUs for parallel processing (default: 1)')
    
    args = parser.parse_args()
    
    # Create metrics calculator and run
    mt = Metrics_scaffold(
        type_cluster=args.type_cluster,
        type_scaffold=args.type_scaffold,
        generator_name=args.generator,
        receptor=args.receptor,
        data_folder=args.data_folder,
        ncpus=args.ncpus
    )
    result = mt.calculate()
    
    if not result.empty:
        print("\n" + "="*60)
        print("FINAL RESULTS:")
        print("="*60)
        print(result.to_string(index=False))


if __name__ == "__main__":
    main()