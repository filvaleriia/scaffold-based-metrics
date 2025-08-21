# Molecular Generator – REINVENT

This repository contains the code and configuration used for running the REINVENT molecular generator, as described in the following publication:

- 📄 **Publication**: [Reinvent 4: Modern AI–driven generative molecule design](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00812-5)  
- 💻 **GitHub (original repository)**: [MolecularAI/REINVENT4](https://github.com/MolecularAI/REINVENT4/tree/main)

---

## 🛠️ Requirements

Before running the generator, make sure to install all required dependencies and clone REINVENT from GitHub.


## 🧪 Workflow Overview

The molecular generation process consists of three main phases:

1. **Transfer Learning (TL)** – adapt a pretrained model to the relevant chemical space.
2. **Staged Learning (SL)** – fine-tune the model to meet specific molecular properties.
3. **Sampling** – generate virtual compounds using the optimized model.

---

## 📥 Input Data

Input molecules are provided for different receptors and split types (similarity or dissimilarity). You can find them in: `data/input_recall_sets/{receptor}/cIS_{receptor}_{type_cluster}_{number}.csv`


- `{receptor}`: Name of the target receptor (e.g., Glucocorticoid_receptor, Leukocyte_elastase)
- `{type_cluster}`: `sim` or `dis`
- `{number}`: Split number from 0 to 4

To convert this input into REINVENT-compatible format, run:

```bash
python data_create_REINVENT.py
```
Processed data will be saved into: `input_set_REINVENT/`

---
## 🔁 Step 1: Transfer Learning

Transfer learning is used to adjust the model to the chemical space of interest. It is configured using a `.toml` file (e.g., `transfer_learning.toml`).

To run transfer learning:

```bash
reinvent transfer_learning.toml
```

- The model is saved after each epoch.
- After training, examine the logs and choose the best-performing epoch for each receptor and cluster combination.
---

## 🧠 Step 2: Staged Learning

Staged learning is used to guide the model toward generating compounds with desired properties.

**Examples of property constraints:**

- Synthetic Accessibility (SA) < 6
- pChEMBL value ≥ 7 (Glucocorticoid receptor)
- pChEMBL value ≥ 6 (Leukocyte elastase)

Use your chosen `staged_learning_train.toml` to run:

```bash
reinvent staged_learning_train.toml
```

Since we want to use our own QSAR model we use in REINVENT external property predictors, pChEMBL values are predicted using a QSAR model (QSPRpred). You can use the script: pchembl_prediction.py

## 🌱 Step 3: Sampling

After staged learning, the **last epoch** is saved as a checkpoint and used to generate a set of virtual compounds.

To generate molecules, use:

```bash
reinvent staged_learning_sampling.toml
```

Final checkpoints are saved here: `checkpoints_SL/{receptor}/`

