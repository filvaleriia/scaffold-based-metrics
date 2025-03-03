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

    cluster_type = 'dis'
    cluster_number = 0
    receptor = 'Glucocorticoid_receptor'


    #qsar
    qsprpred_scorer, sascore = qsar_model(receptor, cluster_type, cluster_number)
    print("qsprpred")
    #creating scoring environment
    environment = creating_scoring_environment(qsprpred_scorer, sascore)

    generated_compaunds(environment, receptor, cluster_number, cluster_type)

    print("Completed_task")