import drugex
import os
import pandas as pd
import numpy as np
import sys

from drugex.data.processing import Standardization
import rdkit
from rdkit import Chem

from drugex.training.monitors import FileMonitor
from drugex.training.scorers.qsprpred import QSPRPredScorer
from drugex.training.scorers.properties import Property
from drugex.training.scorers.modifiers import ClippedScore
from drugex.training.environment import DrugExEnvironment
from drugex.training.rewards import ParetoCrowdingDistance
from drugex.training.scorers.interfaces import Scorer
from drugex.data.corpus.vocabulary import VocSmiles
from drugex.data.datasets import SmilesDataSet
from drugex.training.generators import SequenceRNN
from drugex.training.monitors import FileMonitor
from drugex.data.corpus.corpus import SequenceCorpus
from drugex.data.corpus.vocabulary import  VocSmiles
from drugex.data.datasets import SmilesDataSet
from drugex.data.processing import CorpusEncoder, RandomTrainTestSplitter, Standardization
from utils import smilesToGrid
from drugex.training.explorers import SequenceExplorer
from drugex.training.generators import  SequenceRNN

from qsprpred.models.scikit_learn import SklearnModel


def encoded(receptor, cluster_type, cluster_number):
    #data_set_path = f'data/data_IS/{receptor}'

    df = pd.read_csv(f"cIS_{receptor}_with_p_chembl_{cluster_type}_{cluster_number}.csv", header = None, na_values=('NA', 'nan', 'NaN'))
    df = df.drop([0], axis = 1)
    df.columns = ['SMILES', 'pchembl_value']

    smiles = df.SMILES

    ## Resources
    n_proc = 8 
    batch_size = 512

    # Code
    ## Set paths
    encoded_path = f"data/encoded/rnn/{cluster_type}_{cluster_number}"
    if not os.path.exists(encoded_path):
        os.makedirs(encoded_path)
        

    #MODELS_PR_PATH = 'pretrained_models/DrugEx_v2_PT_Papyrus05.5'
    VOCAB_FILE_PR = 'Papyrus05.5_smiles_rnn_PT.vocab'
    vocabulary = VocSmiles.fromFile(VOCAB_FILE_PR, encode_frags=False)


    ## Standardize SMILES
    standardizer = Standardization(n_proc=n_proc)
    smiles = standardizer.apply(smiles)

    ## Set and run encoder
    encoder = CorpusEncoder(
                SequenceCorpus,
                {
                    'vocabulary': vocabulary, 
                    'update_voc': False, 
                    'throw': True 

                },
                n_proc=n_proc,
                chunk_size=batch_size
    )

    data_collector = SmilesDataSet(f"{encoded_path}/cluster_{cluster_type}_{cluster_number}.tsv", rewrite=True)
    encoder.apply(smiles, collector=data_collector)

    splitter = RandomTrainTestSplitter(0.1, 1e4)
    train, test = splitter(data_collector.getData())
    for data, name in zip([train, test], ['train', 'test']):
        pd.DataFrame(data, columns=data_collector.getColumns()).to_csv(os.path.join(encoded_path, f'cluster_{cluster_type}_{cluster_number}_{name}.tsv'), header=True, index=False, sep='\t')

    vocabulary.toFile(os.path.join(encoded_path, 'pretrained.vocab'))
    print("Encoded end")


def finetuning(receptor, cluster_type, cluster_number):
    print("Start Finetuning")

    ## Pretrained model paths:
    MODEL_FILE_PR = f'Papyrus05.5_smiles_rnn_PT.pkg'
    VOCAB_FILE_PR = f'Papyrus05.5_smiles_rnn_PT.vocab'

    #encoded set path:
    encoded_dir = f'{cluster_type}_{cluster_number}'
    encoded_name = f'cluster_{cluster_type}_{cluster_number}'
    data_set_train_PATH = os.path.join(encoded_dir, f'{encoded_name}_train.tsv')
    data_set_test_PATH = os.path.join(encoded_dir, f'{encoded_name}_test.tsv')

    ## Model parameters
    N_EPOCHS = 1000
    PATIENCE = 100
    BATCH_SIZE = 512

    GPUS = [0,1,2]

    #CODE
    #set path for finutuned
    MODELS_FT_PATH = f'data/finetuned/rnn/{receptor}/{cluster_type}_{cluster_number}'
    if not os.path.exists(MODELS_FT_PATH):
        os.makedirs(MODELS_FT_PATH)

    ## Set pretained model
    vocabulary = VocSmiles.fromFile(VOCAB_FILE_PR, encode_frags=False)

    pretrained = SequenceRNN(voc=vocabulary, is_lstm=True, use_gpus=GPUS)
    pretrained.loadStatesFromFile(MODEL_FILE_PR)

    ## Set encoded test and train sets
    train = SmilesDataSet(data_set_train_PATH, voc=vocabulary)
    train.voc = vocabulary
    train_loader = train.asDataLoader(batch_size=BATCH_SIZE)

    test = SmilesDataSet(data_set_test_PATH, voc=vocabulary)
    test.voc = vocabulary
    valid_loader = test.asDataLoader(batch_size=BATCH_SIZE)

    ## Finetune model
    monitor = FileMonitor(
        f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}', 
        save_smiles=True,
        reset_directory=True,
    )

    pretrained.fit(train_loader, valid_loader, epochs=N_EPOCHS, patience=PATIENCE, monitor=monitor)

    print("Finetuning done.")

    vocabulary.toFile(f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}.vocab')


class ModelScorer(QSPRPredScorer):    
    
    def getScores(self, mols, frags=None):
        """
        Processes molecules and returns a score for each (i.e. a QSAR model prediction).
        If a molecule is invalid, it returns 0 as its score.
        """
        valid_mols = []
        scores = np.zeros(len(mols))  # Default to zero for all

        for i, mol in enumerate(mols):
            if isinstance(mol, str):
                mol = Chem.MolFromSmiles(mol)
            if mol is None:
                continue  # Keep score as 0
            valid_mols.append(mol)

        if valid_mols:
            valid_scores = super().getScores(valid_mols, frags=None)
            scores[: len(valid_scores)] = valid_scores.reshape(len(valid_scores),)

        return scores



def qsar_model(receptor, cluster_type, cluster_number):
    print("QSAR")
    
    predictor = SklearnModel.fromFile(
        os.path.join(f"data/qsprpred/{receptor}/cIS_{receptor}_with_p_chembl_{cluster_type}_{cluster_number}", f"cIS_{receptor}_with_p_chembl_{cluster_type}_{cluster_number}_meta.json")
    )
    
    qsprpred_scorer = ModelScorer(predictor)

    sascore = Property(
        "SA")

    predictor_modifier = ClippedScore(lower_x=4, upper_x=9)
    sascore_modifier = ClippedScore(lower_x=6, upper_x=2)

    qsprpred_scorer.setModifier(predictor_modifier)
    sascore.setModifier(sascore_modifier)

    return qsprpred_scorer, sascore


def creating_scoring_environment(qsprpred_scorer, sascore):
    print("Create env")
    scorers = [
    qsprpred_scorer,
    sascore
    ]
    thresholds = [
        0.6,  #threshold for pchembl to 7
        0.1 
    ]

    environment = DrugExEnvironment(scorers, thresholds, reward_scheme=ParetoCrowdingDistance())

    return environment


def reinforcment(receptor, cluster_type, cluster_number, environment):
    print("Reinforcment")
    # Pretrained model paths:

    MODEL_FILE_PR = f'Papyrus05.5_smiles_rnn_PT.pkg'

    # Finetuned model path:
    MODELS_FT_PATH = f'data/models/finetuned/rnn/{receptor}/{cluster_type}_{cluster_number}'
    MODEL_FILE_FT = f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}.pkg'

    # Vocab path
    VOCAB_FILE = f'Papyrus05.5_smiles_rnn_PT.vocab'

    # Encoded sets path:
    encoded_dir = f'data/encoded/rnn/{cluster_type}_{cluster_number}'
    encoded_name = f'cluster_{cluster_type}_{cluster_number}'
    data_set_train_PATH = os.path.join(encoded_dir, f'{encoded_name}_train.tsv')
    data_set_test_PATH = os.path.join(encoded_dir, f'{encoded_name}_test.tsv')


    ## Model parameters
    N_EPOCHS = 1000
    PATIENCE = 100
    BATCH_SIZE = 512 

    ## Resources
    GPUS = [0,1,2]

    # Code
    # Set paths
    MODELS_RL_PATH = f'data/models/RL/rnn/{receptor}/{cluster_type}_{cluster_number}'
    if not os.path.exists(MODELS_RL_PATH):
        os.makedirs(MODELS_RL_PATH)

    # Set models
    vocabulary = VocSmiles.fromFile(VOCAB_FILE, encode_frags=False)

    pretrained = SequenceRNN(voc=vocabulary, is_lstm=True, use_gpus=GPUS) 
    pretrained.loadStatesFromFile(MODEL_FILE_PR)

    finetuned = SequenceRNN(voc=vocabulary, is_lstm=True, use_gpus=GPUS)
    finetuned.loadStatesFromFile(MODEL_FILE_FT)



    explorer = SequenceExplorer(
            agent=pretrained, 
            env=environment, 
            mutate=finetuned, 
            epsilon=0.1, 
            use_gpus=GPUS,
        )


    print("TADYYY")
    monitor = FileMonitor(f"{MODELS_RL_PATH}/agent", save_smiles=True, reset_directory=True)
    explorer.fit(
        monitor=monitor, 
        epochs=N_EPOCHS,
        patience=PATIENCE
    )
    


def generated_compaunds(environment, receptor, cluster_number, cluster_type):
    print("Generated_compounds")
    # RL path
    MODELS_RL_PATH = f'data/models/RL/rnn/{receptor}/{cluster_type}_{cluster_number}'

    # Pretrained model paths:
    MODEL_FILE_PR = f'Papyrus05.5_smiles_rnn_PT.pkg'

    # Vocab path
    VOCAB_FILE_PR = f'Papyrus05.5_smiles_rnn_PT.vocab'


    ## Resources
    GPUS = [0,1,2]


    vocabulary = VocSmiles.fromFile(VOCAB_FILE_PR, encode_frags=False)


    reinforced = SequenceRNN(voc=vocabulary, is_lstm=True, use_gpus=GPUS)
    reinforced.loadStatesFromFile(f'{MODELS_RL_PATH}/agent.pkg')

    generated = pd.DataFrame()

    generated = reinforced.generate(num_samples=200000, batch_size=1024, n_proc=8, evaluator=environment, raw_scores=True)

    print("Generated compounds done")
    generated.to_csv(f'cOS_DrugEx_epsilon_0.6_{cluster_type}_{cluster_number}_all_columns.csv', index = False)
    print("CSV savings")

    return generated



if '__main__':

    cluster_type = 'dis'
    cluster_number = 0
    receptor = 'Glucocorticoid_receptor'
    #save encoded file to folder 
    encoded(receptor, cluster_type, cluster_number)

    #finetuning
    finetuning(receptor, cluster_type, cluster_number)

    #qsar
    qsprpred_scorer, sascore = qsar_model(receptor, cluster_type, cluster_number)

    #creating scoring environment
    environment = creating_scoring_environment(qsprpred_scorer, sascore)

    #reinforcment/training of model
    reinforcment(receptor, cluster_type, cluster_number, environment)
  

    generated_compaunds(environment, receptor, cluster_number, cluster_type)


