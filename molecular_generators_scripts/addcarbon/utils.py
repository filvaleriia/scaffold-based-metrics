import uuid
from functools import partial
from multiprocessing import Pool
from time import gmtime, strftime

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import roc_auc_score


def timestamp(adduuid=False):
    s = strftime("%Y-%m-%d_%H:%M:%S", gmtime())
    if adduuid:
        s = s + '_' + uuid.uuid4().hex
    return s



