import pandas as pd
import os


receptor = 'Glucocorticoid_receptor'
for type in ['dis', 'sim']:
    for number in range(5):
        df = pd.read_csv(f"../../data/input_recall_sets/{receptor}/cIS_{receptor}_with_p_chembl_{type}_{number}.csv", header = None)
        df.columns = ['id','smiles','p_chembl']
        df = df.drop(columns= ['id'])
        data = df.sample(frac=1)
        n_head = len(data) // 5
        n_tail = len(df) - n_head
               
        
        train, validation = data.head(n_tail), data.tail(n_head)
        
        folder = f"input_set_REINVENT/{receptor}"
        if not os.path.exists(folder):
            os.makedirs(folder)
        train.to_csv(f"{folder}/cIS_GR_{type}_{number}.smi", sep="\t", index=False, header=False)
        validation.to_csv(f"{folder}/cIS_GR_{type}_{number}_valid.smi", sep="\t", index=False, header=False)