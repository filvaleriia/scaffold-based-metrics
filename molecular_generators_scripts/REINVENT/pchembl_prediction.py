from qsprpred.models import SklearnModel
import sys
import json
from rdkit import Chem
import numpy as np


number = 1
type = 'sim'

def predict_pchembl(smiles):
    model = SklearnModel.fromFile(f"qsprpred/Glucocorticoid_receptor/cIS_Glucocorticoid_receptor_with_p_chembl_{type}_{number}/cIS_Glucocorticoid_receptor_with_p_chembl_{type}_{number}_meta.json")
    return model.predictMols(smiles)


smilies = [smiles.strip() for smiles in sys.stdin]

scores = predict_pchembl(smilies)

dict_ = {"version": 1, "payload": {"predictions": np.reshape(scores,len(smilies)).tolist()}}
print(json.dumps(dict_))


