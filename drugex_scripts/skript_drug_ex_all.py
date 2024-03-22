#zatim jenom skript pro generovani sloucenin na jiz netranovanem generatoru

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

from qsprpred.models import SklearnModel


def encoded(receptor, cluster_type, cluster_number):
    #path to data sets
    data_set_path = f'recall_metrics/data/{receptor}'


    ## Fragmentation method:
    frag_method = 'brics'

    n_proc = 8
    ## Retrieve dataset
    df = pd.read_csv(f"{data_set_path}/IS_{receptor}_with_p_chembl_{cluster_type}_{cluster_number}.csv", header = None, na_values=('NA', 'nan', 'NaN'))
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
    encoded_dir = f'recall_metrics/data/{receptor}/encoded/graph'
    encoded_name = f'cluster_{cluster_type}_{cluster_number}'

    ## Model parameters
    N_EPOCHS = 1000
    PATIENCE = 100
    BATCH_SIZE = 512

    GPUS = [0,1,2]

    #set path for finutuned
    MODELS_FT_PATH = f'recall_metrics/models/finetuned/graph/{receptor}/{cluster_type}_{cluster_number}'
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
        name=f'Glucocorticoid_{cluster_type}_{cluster_number}_RandomForestClassifier',
        base_dir=f"recall_metrics/data/{receptor}/qsar"
        )

    qsprpred_scorer = QSPRPredScorer(
        predictor)

    sascore = Property(
        "SA")

    predictor_modifier = ClippedScore(lower_x=0.2, upper_x=0.8)
    sascore_modifier = ClippedScore(lower_x=5, upper_x=2)

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
    MODELS_FT_PATH = f'recall_metrics/models/finetuned/graph/{receptor}/{cluster_type}_{cluster_number}'

    # Vocab path
    VOCAB_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT.vocab'
    VOCAB_FILE_FT = f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}.vocab'

    # Encoded sets path:
    encoded_dir = f'recall_metrics/data/{receptor}/encoded/graph'
    encoded_name = f'cluster_{cluster_type}_{cluster_number}'

    ## Model parameters
    N_EPOCHS = 1000
    PATIENCE = 100
    BATCH_SIZE = 512

    ## Resources
    GPUS = [0,1,2]

    # Code
    # Set paths
    MODELS_RL_PATH = f'recall_metrics/models/RL/graph/{receptor}/{cluster_type}_{cluster_number}'
    if not os.path.exists(MODELS_RL_PATH):
        os.makedirs(MODELS_RL_PATH)

    # Set models
    vocabulary_PR = VocGraph.fromFile(VOCAB_FILE_PR)
    vocabulary_FT = VocGraph.fromFile(VOCAB_FILE_FT)

    pretrained = GraphTransformer(voc_trg=vocabulary_PR, use_gpus=GPUS)
    pretrained.loadStatesFromFile(MODEL_FILE_PR)

    finetuned = GraphTransformer(voc_trg=vocabulary_FT, use_gpus=GPUS)
    finetuned.loadStatesFromFile(f'{MODELS_FT_PATH}/cluster_{cluster_type}_{cluster_number}.pkg')

    explorer = FragGraphExplorer(agent=pretrained, env=environment, mutate=finetuned, epsilon=0.2, use_gpus=GPUS)

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



def generated_compaunds(environment, receptor, cluster_number, cluster_type):
    print("Generated_compounds")
    # RL path
    MODELS_RL_PATH = f'recall_metrics/models/RL/graph/{receptor}/{cluster_type}_{cluster_number}'

    # Encoded sets path:
    encoded_dir = f'recall_metrics/data/{receptor}/encoded/graph'
    encoded_name = f'cluster_{cluster_type}_{cluster_number}'

    # Pretrained model paths:
    MODELS_PR_PATH = 'Papyrus05.5_graph_trans_PT'
    MODEL_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT.pkg'

    # Vocab path
    VOCAB_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT.vocab'


    ## Resources
    GPUS = [0,1,2]

    train = GraphFragDataSet(f"{encoded_dir}/{encoded_name}_train.tsv")
    test = GraphFragDataSet(f"{encoded_dir}/{encoded_name}_test.tsv")

    vocabulary = VocGraph.fromFile(VOCAB_FILE_PR)
    reinforced = GraphTransformer(voc_trg=vocabulary, use_gpus=GPUS)
    reinforced.loadStatesFromFile(f'{MODELS_RL_PATH}/agent.pkg')

    generated = pd.DataFrame()
    #generated = reinforced.generate(input_dataset=train_from_file, num_samples=10, evaluator=environment, raw_scores=True)
    generated = reinforced.generate(input_dataset=train, num_samples=200000, evaluator=environment, raw_scores=True)
    #generated = reinforced.generate(input_dataset=train_from_file, num_samples=20, evaluator=environment, raw_scores=True)
    print("Generated compounds done")
    generated.to_pickle(f'cOS_DrugEx_{cluster_type}_{cluster_number}_all_columns.pkl')
    print("Pickle savings")
    generated.to_csv(f'cOS_DrugEx_{cluster_type}_{cluster_number}_all_columns.csv', index = False)
    print("CSV savings")
    return generated

if '__main__':

    cluster_type = 'sim'
    cluster_number = 1
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
    
    returr = generated_compaunds(environment, receptor, cluster_number, cluster_type)

    print("Completed_task")
