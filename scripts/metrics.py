'''This script contains the count for three important matrics: TUPOR, SASY, ASER'''


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
    def __init__(self, type_cluster, scaffold_type, generator_name, number_of_calculation, receptor, save_options):
        self.output_set = []
        self.recall_set = []
        self.output_set_scaffolds = []
        self.recall_set_scaffolds = []
        self.type_cluster = type_cluster
        self.unique_output_set = []
        self.unique_recall_set = []
        self.count_metrics = []
        self.scaffold_type = scaffold_type
        self.generator_name = generator_name
        self.number_of_calculation = number_of_calculation
        self.save_options = save_options
        self.receptor = receptor
        self.results = pd.DataFrame()
    

    def load(self,filepath_output_set, filepath_recall_set):
        '''Load tha main dataset needed for calculation matrics, the first dataset is the Output Set and the second one is the Recall Set. 
        All datasets add to constructor self.output_set and self.recall_set'''

        output_set = []
        recall_set = []
        print(filepath_output_set)
        with open(filepath_output_set, 'r') as f:
            for line in f.readlines():
                output_set.append(line)
        self.output_set = pd.DataFrame(output_set)


        with open(filepath_recall_set, 'r') as f:
            for line in f.readlines():
                recall_set.append(line)
        self.recall_set = pd.DataFrame(recall_set)


    def scaffolds_unique(self,output_set,recall_set):
        '''This function only return unique scaffolds in Output Set and Recall Set'''
        return output_set[0].unique(), recall_set[0].unique()
    

    def main_function_return(self,output_set, recall_set):
        '''The main function for calculation all metrics'''

        '''Convert all metrics to scaffold. Possible scaffolds: csk and murcko'''
        self.recall_set_scaffolds = convert_to_scaffold_(recall_set,self.scaffold_type)
        self.output_set_scaffolds = convert_to_scaffold_(output_set, self.scaffold_type)

        '''Create unique dataset and save to constructor like self.unique_output_set and self.unique_recall_set'''
        self.unique_output_set, self.unique_recall_set = self.scaffolds_unique(self.output_set_scaffolds, self.recall_set_scaffolds)

        '''Calculate the occurance of scaffolds'''
        df = add_columns_same_like_input_function(self.output_set_scaffolds, self.unique_recall_set)

        '''Calculate the individual metrics'''
        ns = len(self.unique_output_set)
        ss = len(self.output_set_scaffolds)
        sesy = ns/ss
        aser = df[1].sum()/ss

        self.count_metrics = df.copy()
        self.count_metrics.columns = ['unique_scaffold_recall','count_of_occurance', 'uniq_occurance']

        '''If generator doesn't generate unique active compounds the number of TUPOR'''
        try:
            tupor = df[2].value_counts()[1]/len(df)
            tupor_text = f"{df[2].value_counts()[1]}/{len(df)}"
        except:
            tupor = 0
            tupor_text = f"{0}/{len(df)}"


        tRS = df[1].sum()

        '''Return individual metrics and next add to pandas data frame'''
        return self.type_cluster , ns,ss,tupor_text,tupor,sesy,aser,tRS


    def save_function(self):
        '''Save function to folder 'data/results'. If the folder doesn't exist -> the folder will be created'''
        print("Save")
        path_to_folder = f'data/results/{self.receptor}'
        
        '''Check if folder is existing'''
        if not os.path.exists(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}"):
            os.makedirs(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}")
        
        self.count_metrics.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/count_of_occurance_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)

        self.results.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/metrics_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)

        self.output_set_scaffolds.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/scaffolds_of_output_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)

        self.recall_set_scaffolds.to_csv(f"{path_to_folder}/{self.scaffold_type}_scaffolds/{self.type_cluster}/{self.generator_name}/scaffolds_of_recall_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)


    def calculate_metrics_with_return(self):
        '''Calculation all 7 metrics'''

        res = self.main_function_return(self.output_set, self.recall_set)

        results = pd.DataFrame(columns = ['type_cluster','uniq_scaffolds','set_size','TUPOR_','TUPOR',\
                                 'SESY','ASER', 'tRS'])
        results.loc[len(results)] = res
        results.insert(loc=0, column='name', value=[f"{self.generator_name}_{self.number_of_calculation}"])
        results.insert(loc=2, column='scaffold', value=[self.scaffold_type])

        self.results = results

        if self.save_options == True:
            self.save_function()
        return results
    

    


    
       



















