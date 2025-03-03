import drugex
import os
import pandas as pd
import numpy as np

from drugex.data.processing import Standardization
from drugex.data.fragments import FragmentCorpusEncoder
from drugex.data.fragments import GraphFragmentEncoder, FragmentPairsSplitter
from drugex.molecules.converters.fragmenters import Fragmenter
from drugex.data.corpus.vocabulary import VocGraph
from drugex.data.datasets import GraphFragDataSet
from drugex.data.corpus.vocabulary import VocGraph
from drugex.training.generators import GraphTransformer
from drugex.training.monitors import FileMonitor
from drugex.training.scorers.qsprpred import QSPRPredScorer
from drugex.training.scorers.properties import Property
from drugex.training.scorers.modifiers import ClippedScore
from drugex.training.environment import DrugExEnvironment
from drugex.training.rewards import ParetoCrowdingDistance
from drugex.training.explorers import FragGraphExplorer

from utils import smilesToGrid

from qsprpred.models.scikit_learn import SklearnModel


def encoded(receptor, cluster_type, cluster_number):
    #path to data sets
    data_set_path = f'data/data_IS/{receptor}'


    ## Fragmentation method:
    frag_method = 'brics'

    n_proc = 8
    ## Retrieve dataset
    df = pd.read_csv(f"{data_set_path}/cIS_{receptor}_with_p_chembl_{cluster_type}_{cluster_number}.csv", header = None, na_values=('NA', 'nan', 'NaN'))
    df = df.drop([0], axis = 1)
    df.columns = ['SMILES', 'pchembl_value']

    smiles = df.SMILES

    ## Standardize SMILES
    standardizer = Standardization(n_proc=n_proc)
    smiles = standardizer.apply(smiles)
    smiles = set(smiles)
    
    ## Set and run encoder
    encoder = FragmentCorpusEncoder(
        fragmenter = Fragmenter(4, 4, frag_method), 
        encoder = GraphFragmentEncoder(VocGraph(n_frags=4)),
        pairs_splitter = FragmentPairsSplitter(), 
        n_proc=n_proc,
    )

    #set a path
    encoded_path = f"{data_set_path}/encoded/graph"
    if not os.path.exists(encoded_path):
        os.makedirs(encoded_path)
        
    # create empty data sets (we have to specify a path to a file where the data set will be saved)
    train = GraphFragDataSet(f"{encoded_path}/cluster_{cluster_type}_{cluster_number}_train.tsv", rewrite=True)
    test = GraphFragDataSet(f"{encoded_path}/cluster_{cluster_type}_{cluster_number}_test.tsv", rewrite=True)

    encoder.apply(list(smiles), encodingCollectors=[test, train])
    print("Encoded end")


def finetuning(receptor, cluster_type, cluster_number):
    print("Start Finetuning")
    ## Pretrained model paths:
    MODELS_PR_PATH = 'Papyrus05.5_graph_trans_PT'
    MODEL_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT.pkg'
    VOCAB_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT.vocab'

    #encoded set path:
    encoded_dir = f'data/data_IS/{receptor}/encoded/graph'
    encoded_name = f'cluster_{cluster_type}_{cluster_number}'

    ## Model parameters
    N_EPOCHS = 1000
    PATIENCE = 100 
    BATCH_SIZE = 512

    GPUS = [0,1,2]

    #set path for finutuned
    MODELS_FT_PATH = f'data/models/finetuned/graph/{receptor}/{cluster_type}_{cluster_number}'
    if not os.path.exists(MODELS_FT_PATH):
        os.makedirs(MODELS_FT_PATH)

    ## Set pretained model
    vocabulary = VocGraph.fromFile(VOCAB_FILE_PR)

    pretrained = GraphTransformer(voc_trg=vocabulary, use_gpus=GPUS) 
    pretrained.loadStatesFromFile(MODEL_FILE_PR)

    ## Set encoded test and train sets
    train = GraphFragDataSet(f"{encoded_dir}/{encoded_name}_train.tsv")
    test = GraphFragDataSet(f"{encoded_dir}/{encoded_name}_test.tsv")

    ## Finetune model
    monitor = FileMonitor(
        f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}', 
        save_smiles=True,
        reset_directory=True,
    )

    pretrained.fit(
        train.asDataLoader(BATCH_SIZE),
        test.asDataLoader(BATCH_SIZE),
        epochs=N_EPOCHS,
        patience=PATIENCE, 
        monitor=monitor 
    )

    print("Finetuning done.")

    vocabulary.toFile(f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}.vocab')


def qsar_model(receptor, cluster_type, cluster_number):
    print("QSAR")
    predictor = SklearnModel(
        name=f'cIS_Glucocorticoid_receptor_with_p_chembl_{cluster_type}_{cluster_number}',
        base_dir=f"data/qsprpred/{receptor}"
        )

    qsprpred_scorer = QSPRPredScorer(
        predictor)

    sascore = Property(
        "SA")

    predictor_modifier = ClippedScore(lower_x=0.2, upper_x=0.8)
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
        0.5, 
        0.1 
    ]

    environment = DrugExEnvironment(scorers, thresholds, reward_scheme=ParetoCrowdingDistance())

    return environment


def reinforcment(receptor, cluster_type, cluster_number, environment):
    print("Reinforcment")
    # Pretrained model paths:
    MODELS_PR_PATH = 'Papyrus05.5_graph_trans_PT'
    MODEL_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT.pkg'

    # Finetuned model path:
    MODELS_FT_PATH = f'data/models/finetuned/graph/{receptor}/{cluster_type}_{cluster_number}'

    # Vocab path
    VOCAB_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT.vocab'
    VOCAB_FILE_FT = f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}.vocab'

    # Encoded sets path:
    encoded_dir = f'data/data_IS/{receptor}/encoded/graph'
    encoded_name = f'cluster_{cluster_type}_{cluster_number}'

    ## Model parameters
    N_EPOCHS = 1000
    PATIENCE = 100 
    BATCH_SIZE = 512 

    ## Resources
    GPUS = [0,1,2]

    # Code
    # Set paths
    MODELS_RL_PATH = f'data/models/RL/graph/{receptor}/{cluster_type}_{cluster_number}'
    if not os.path.exists(MODELS_RL_PATH):
        os.makedirs(MODELS_RL_PATH)

    # Set models
    vocabulary_PR = VocGraph.fromFile(VOCAB_FILE_PR)
    vocabulary_FT = VocGraph.fromFile(VOCAB_FILE_FT)

    pretrained = GraphTransformer(voc_trg=vocabulary_PR, use_gpus=GPUS)
    pretrained.loadStatesFromFile(MODEL_FILE_PR)

    finetuned = GraphTransformer(voc_trg=vocabulary_FT, use_gpus=GPUS)
    finetuned.loadStatesFromFile(f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}.pkg')

    explorer = FragGraphExplorer(agent=finetuned, env=environment, mutate=pretrained, epsilon=0.2, use_gpus=GPUS)

    train = GraphFragDataSet(f"{encoded_dir}/{encoded_name}_train.tsv")
    test = GraphFragDataSet(f"{encoded_dir}/{encoded_name}_test.tsv")

    print("TADYYY")
    monitor = FileMonitor(f"{MODELS_RL_PATH}/agent", save_smiles=True, reset_directory=True)
    explorer.fit(
        train.asDataLoader(batch_size=BATCH_SIZE), 
        test.asDataLoader(batch_size=BATCH_SIZE), 
        monitor=monitor, 
        epochs=N_EPOCHS,
        patience=PATIENCE)





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
    print("qsprpred")
    #creating scoring environment
    environment = creating_scoring_environment(qsprpred_scorer, sascore)

    #reinforcment/training of model
    reinforcment(receptor, cluster_type, cluster_number, environment)
    print("reinforcment")
    


    print("Completed_task")
