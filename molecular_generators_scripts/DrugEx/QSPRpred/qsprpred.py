# Imports
import os
import pandas as pd
from qsprpred.data.descriptors.fingerprints import MorganFP
from qsprpred.data.descriptors.sets import DrugExPhyschem
from qsprpred.data.sampling.splits import ScaffoldSplit
from qsprpred.data import QSPRDataset
from qsprpred.models import OptunaOptimization, TestSetAssessor, CrossValAssessor, FileMonitor
from qsprpred.models.scikit_learn import SklearnModel
import xgboost as xgb

# Paths
work_dir = os.getcwd()
results_dir = os.path.join(work_dir, 'results/')
os.makedirs(results_dir, exist_ok=True)

final_models = os.path.join(work_dir, 'final_models/')
os.makedirs(final_models, exist_ok=True)

data_path = os.path.join(work_dir, 'data_Leukocyte_elastase/')
print(data_path)

input_sets = [
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_dis_0.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_dis_1.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_dis_2.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_dis_3.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_dis_4.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_sim_0.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_sim_1.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_sim_2.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_sim_3.csv'),
    os.path.join(data_path, 'cIS_Leukocyte_elastase_with_p_chembl_sim_4.csv'),
]

 

print(input_sets)
inactive_set = [
    os.path.join(data_path, 'Leukocyte_elastase_non_active_compounds.csv'),
]

df_results = pd.DataFrame(columns=['name', 'r2', 'rmse'])
df_results_list = [df_results]

for input_set in input_sets:
    print(input_set)
    # Prepare datasets
    df_active = pd.read_csv(input_set, index_col=0, header=None)
    df_active.columns = ['canonical_smiles', 'pchembl_value']
    df_inactive = pd.read_csv(inactive_set[0])
    df = pd.concat([df_active, df_inactive])

    dataset = QSPRDataset(
        name = os.path.splitext(os.path.basename(input_set))[0],
        store_dir = work_dir,
        df = df,
        target_props = [{
            "name": "pchembl_value", 
            "task": "REGRESSION",
        }],
        smiles_col = "canonical_smiles",
        random_state = 42,
        n_jobs = 16,
    )

    descriptors = [
        MorganFP(radius=3, nBits=2048), 
        DrugExPhyschem(),
    ]

    dataset.prepareDataset(
        split=ScaffoldSplit(test_fraction=0.2, dataset=dataset),
        feature_calculators=descriptors,
        recalculate_features=True,
    )

    dataset.save()

    # Modelling
    fixed_paramaters = {
        "objective": "reg:squarederror",
        "verbosity": 0,
        "n_jobs": 16
    }

    search_space = {
        "n_estimators": ['int', 100, 2500],
        "learning_rate": ['float', 1e-3, 0.1],
        "max_depth": ['int', 1, 10],
        "subsample": ['float', 0.05, 1.0],
        "colsample_bytree": ['float', 0.05, 1.0],
        "min_child_weight": ['int', 1, 20],
    }

    model = SklearnModel(
        name=os.path.splitext(os.path.basename(input_set))[0],
        base_dir=os.path.join(work_dir, 'models'),
        alg=xgb.XGBRegressor, 
        parameters=fixed_paramaters,
    )

    gridsearcher = OptunaOptimization(
        n_trials=300,
        param_grid=search_space,
        model_assessor=CrossValAssessor(scoring='r2', split=ScaffoldSplit(n_folds=5)),
        monitor=FileMonitor()
    )

    gridsearcher.optimize(model, dataset)

    df_results_list.append(pd.DataFrame([{
            'name': os.path.splitext(os.path.basename(input_set))[0], 
            'r2': TestSetAssessor('r2', use_proba=False)(model, dataset)[0],
            'rmse': -1 * TestSetAssessor('neg_root_mean_squared_error', use_proba=False)(model, dataset)[0],
        }
    ]))

    _ = model.fitDataset(dataset, monitor=FileMonitor())

    model.save(os.path.join(work_dir, 'final_models/', os.path.splitext(os.path.basename(input_set))[0]))
    print(model.save())

df_results = pd.concat(df_results_list, ignore_index=True)
df_results.to_csv(os.path.join(work_dir, 'results/model_results.csv'), index=False)
print(df_results)