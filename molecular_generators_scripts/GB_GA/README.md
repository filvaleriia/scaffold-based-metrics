# Molecular Generator GB_GA

This repository contains the code for the genetic algorithm-based molecular generator described in the following publication:

- 📄 **Publication**: [Jensen et al., RSC Chemical Science (2019)](https://pubs.rsc.org/en/content/articlelanding/2019/sc/c8sc05372c)  
- 💻 **GitHub (original repository)**: [https://github.com/jensengroup/GB_GA](https://github.com/jensengroup/GB_GA/tree/master)

---

## 🛠️ Requirements

Before running the generator, please install the required packages. You will need:
- `rdkit`
- `numpy`
- `pandas`

You can install all dependencies using:

```bash
pip install -r requirements.txt
```

---

## ⚙️ Parameters and Customization

You can modify several parameters directly in the script to control the behavior of the generator, such as:

- **Mutation rate**: 
- **Population size**
- **Number of generations**
- **Mating pool size**
- **Scoring function**
- **Size of molecules**

To adjust them, open `GB_GA.py` and modify the corresponding variables in the configuration section.

---

## 📂 Input Data

Each run requires an input set of molecules. These are located under: 
`data/input_recall_sets/{receptor}/`

Each file follows the naming pattern:
`cIS_{receptor}{type_cluster}{number}.csv`


Where:
- `receptor` is the target (e.g., `Glucocorticoid_receptor`, `Leukocyte_elastase`)
- `type_cluster` is either `dis` or `sim`
- `number` is in the range `0–4`

Make sure to provide the correct path when running the script.

---

## 🧪 Example Run

Generate molecules for Glucocorticoid receptor(name of receptor you can change inside script) using dissimilarity split 0:

```bash
python3 GB_GA.py dis 0
```



