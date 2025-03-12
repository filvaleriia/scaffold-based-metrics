# 🧪 Molecular Generator – AddCarbon

The AddCarbon molecular generator introduces random modifications by replacing atoms with carbon atoms in the input set of molecules.

- 📄 **Paper**: [Why Do These Molecules Look So Similar?](https://www.sciencedirect.com/science/article/pii/S1740674920300159)  
- 💻 **GitHub**: [ml-jku/mgenerators-failure-modes – addcarbon.py](https://github.com/ml-jku/mgenerators-failure-modes/blob/master/addcarbon.py)

---

## ⚙️ Installation

To run this generator, you need the following packages:  
`rdkit`, `numpy`, `sklearn` and `python-Levenshtein`.

You can either install them into your existing environment or use the `requirements.txt`:

```bash
pip install -r requirements.txt
```

## 🚀 Usage

To run the script, simply adjust the input file name inside `addcarbon.py`, and then execute:

```bash
python3 addcarbon.py
```

You can use the input file located at `data/input_recall_sets/{receptor}/cIS_{receptor}_{type_cluster}_{number}.csv` directly without any modifications.
