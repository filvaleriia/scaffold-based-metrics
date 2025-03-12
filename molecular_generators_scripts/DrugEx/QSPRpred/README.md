# QSPRpred

For the DrugEx generator and REINVENT, we created QSAR models using the QSPRpred library. We used the XGBoost algorithm and Optuna optimization to find the best hyperparameters.

- 💻 **GitHub**: [QSPRpred Repository](https://github.com/CDDLeiden/QSPRpred/tree/b981855556fcf91ed1e45895a11f2f2f574366fd)
- 📚 **Documentation**: [QSPRpred Documentation](https://cddleiden.github.io/QSPRpred/docs/)

---

## 🛠️ Running QSPRpred

To run QSPRpred, follow these steps:

1. Install the required dependencies. I used the same environment as for running the DrugEx scripts.

```bash
pip install -r requirements.txt
```

2. Specify the input data in the script.
3. Run the Python script to train the model.

```
python3 qsprpred.py
```

The output models will be saved in the final_models folder.
