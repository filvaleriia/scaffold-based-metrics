# Molecular Generator Molpher

- 📄 **Publication**: [Molecular Generation with Molpher](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-6-7)
- 💻 **GitHub**: [Molpher Repository](https://github.com/lich-uct/molpher-lib)

---

## 🛠️ Installation

Before running the generator, make sure to install all required dependencies and clone Molpher from GitHub.

## 📥 Input Data

Molpher requires specific input data, which includes a starting molecule and a target molecule. The generator will attempt to create a path between these molecules.

The input sets for Molpher can be found in: `data/input_recall_sets/{receptor}/cIS_Molpher_{receptor}_{type_cluster}_{number}.csv`


Where:
- `{receptor}`: Name of the receptor (e.g., Glucocorticoid_receptor, Leukocyte_elastase)
- `{type_cluster}`: The type of cluster (e.g., `sim` for similarity or `dis` for dissimilarity)
- `{number}`: The split number (e.g., `0`, `1`, `2`, etc.)

These files are used as input for Molpher, where each pair consists of a starting molecule and a target molecule, and the model will attempt to create a synthetic path between them based on the given data.


## 🖥️ Running Molpher

To run Molpher, use the following command:

```bash
python3 run_class.py id_pair start_id stop_id start_smiles stop_smiles
```

Where:
- `id_pair`: The identifier for the molecule pair.
- `start_id`: The ID of the starting molecule.
- `stop_id`: The ID of the target molecule.
- `start_smiles`: The SMILES string for the starting molecule.
- `stop_smiles`: The SMILES string for the target molecule.


This command will initiate the path generation between the given starting molecule and the target molecule based on the provided data and parameters.



