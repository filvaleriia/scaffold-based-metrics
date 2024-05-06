'''Import necessary libraries'''

import pandas as pd
import numpy as np
import lzma
import csv
import sys
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric
from rdkit.Chem import Lipinski
from rdkit.Chem.Lipinski import NumAromaticHeterocycles
from rdkit.Chem.Lipinski import NumAliphaticRings
from rdkit.Chem.Lipinski import NumAromaticHeterocycles
from rdkit.Chem.Lipinski import NumAromaticRings
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import Draw

import seaborn as sns
import matplotlib.pyplot as plt





def convert_to_scaffold_(df, scaffold_type):
    '''Function wich is convert compounds to scaffold, two options of scaffold is SCK scaffold and Murcko scaffold'''
    a = []
    print("Convert_")
    if scaffold_type == 'csk':
        for x in range(len(df)):
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(df.loc[x][0]))
            
            try:
                a.append(MurckoScaffoldSmiles(\
                                         Chem.MolToSmiles(MakeScaffoldGeneric\
                                       (Chem.MolFromSmiles(smiles)))))
            except:
                print("Faild to create scaffold_csk")
                print("Index",x)
                print(df.loc[x][0])
                print(smiles)
    elif scaffold_type == 'murcko':
        for x in range(len(df)):
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(df.loc[x][0]))
            try:
                a.append(MurckoScaffoldSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smiles))))
                
            except:
                print("Faild to create scaffold_csk")
                print("Index",x)
                print(smiles)

    dff = pd.DataFrame(data = a)
    return dff


def add_columns_same_like_input_function(df_generated, test_set):
    print("add_columns_same_like_input_function")
    '''Create table with three columns where the first one is 'The uniq compounds of test set', the second one is tRS => 'Number of compounds wich have active scaffold',
    the last one NAS => 'The number of unique active scaffold' '''
    df_generated = pd.DataFrame(df_generated)   
    df = pd.DataFrame()
    df[0] = test_set
    df[1] = int(0)
    for x in range(len(test_set)):
        y = df[0][x]
        df[1][x] = [df_generated[0].value_counts()[y] if y in df_generated[0].unique() else 0][0]       
    
    df[2] = df[1].apply(lambda x : 1 if x > 0 else 0)

    return df


class Metrics:
    def __init__(self,type_cluster,scaffold_type, generator_name, number_of_calculation, save_options):
        self.generated_compounds = []
        self.test_set = []
        self.generated_compounds_scaffolds = []
        self.test_set_scaffolds = []
        self.type_cluster = type_cluster
        self.unique_generated_compounds = []
        self.unique_test_set = []
        self.count_metrics = []
        self.count_metrics_repead = []
        self.scaffold_type = scaffold_type
        self.generator_name = generator_name
        self.number_of_calculation = number_of_calculation
        self.save_options = save_options
        self.results = pd.DataFrame()
    

    def load(self,filepath_generated_com, filepath_test_set):
        '''Load tha main dataset needed for calculation matrics, the first dataset is generated compounds and the second one is the test set. All datasets add to constructor self.generated_compounds and self.test_set'''

        generated_compounds = []
        test_set = []
    
        with open(filepath_generated_com, 'r') as f:
            for line in f.readlines():
                generated_compounds.append(line)
        self.generated_compounds = pd.DataFrame(generated_compounds)


        with open(filepath_test_set, 'r') as f:
            for line in f.readlines():
                test_set.append(line)
        self.test_set = pd.DataFrame(test_set)


    def scaffolds_unique(self,df_generated,df_test_set):
        '''This function only return unique scaffolds in generated set and test set'''
        return df_generated[0].unique(), df_test_set[0].unique()
    

    def main_function_return(self,df_generated, test_set):
        '''The main function for calculation all metrics'''

        '''Convert all metrics to scaffold. Possible scaffolds: csk and murcko'''
        self.test_set_scaffolds = convert_to_scaffold_(test_set,self.scaffold_type)
        self.generated_compounds_scaffolds = convert_to_scaffold_(df_generated, self.scaffold_type)

        '''Create unique dataset and save to constructor like self.unique_generated_compounds and self.unique_test_set'''
        self.unique_generated_compounds, self.unique_test_set = self.scaffolds_unique(self.generated_compounds_scaffolds, self.test_set_scaffolds)

        '''Calculate the occurance of scaffolds'''
        df = add_columns_same_like_input_function(self.generated_compounds_scaffolds, self.unique_test_set)

        df1 = add_columns_same_like_input_function(self.generated_compounds_scaffolds, self.test_set_scaffolds)

        '''Calculate the individual metrics'''
        ns = len(self.unique_generated_compounds)
        ss = len(self.generated_compounds_scaffolds)
        sescy = ns/ss
        sescr = df[2].sum()/ns
        sescry = sescy*sescr
        asescr = df[1].sum()/ss

        self.count_metrics = df.copy()
        self.count_metrics_repead = df1.copy()
        self.count_metrics.columns = ['scaffold_test_tp','count_of_occurance', 'uniq_occurance']
        self.count_metrics_repead.columns = ['scaffold_test_tp_with_repeat','count_of_occurance', 'uniq_occurance']
        
        '''If generator doesn't generate unique active compounds the number of TPRA and TPRAR is 0'''
        try:
            tpra = df[2].value_counts()[1]/len(df)
            tpra_text = f"{df[2].value_counts()[1]}/{len(df)}"
        except:
            tpra = 0
            tpra_text = f"{0}/{len(df)}"

        try:
            tprar = df1[2].value_counts()[1]/len(df1)
            tprar_text = f"{df1[2].value_counts()[1]}/{len(df1)}"
        except:
            tprar = 0
            tprar_text = f"{0}/{len(df1)}"

        tRS = df[1].sum()

        '''Return individual metrics and next add to pandas data frame'''
        return self.type_cluster , ns,ss,tpra_text,tpra, tprar_text,tprar,sescy,sescr, sescry,asescr,tRS


    def save_function(self):
        '''Save function to folder 'data/results'. If the folder doesn't exist -> the folder will be created'''
        print("Save")
        path_to_folder = 'data/results'
        
        '''Check if folder is existing'''
        if not os.path.exists(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}"):
            os.makedirs(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}")
        
        self.count_metrics.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/count_of_occurance_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)

        self.count_metrics_repead.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/count_of_occurance_with_repead_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)

        self.results.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/metrics_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)

        self.generated_compounds_scaffolds.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/scaffolds_of_generated_moleculs_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)

        self.test_set_scaffolds.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/scaffolds_of_test_moleculs_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)


    def calculate_metrics_with_return(self):
        '''Calculation all 7 metrics'''
        res = self.main_function_return(self.generated_compounds, self.test_set)

        results = pd.DataFrame(columns = ['type_cluster','uniq_scaffolds','set_size','tpra_','tpra','tprar_','tprar',\
                                 'sescy','sescr','sescry','asescr', 'tRS'])
        results.loc[len(results)] = res
        results.insert(loc=0, column='name', value=[f"{self.generator_name}_{self.number_of_calculation}"])
        results.insert(loc=2, column='scaffold', value=[self.scaffold_type])

        self.results = results

        if self.save_options == True:
            self.save_function()
        return results
    

    


    
       



















