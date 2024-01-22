'''import important library'''
import dataframe_image as dfi
from pandas.plotting import table
from html2image import Html2Image

from matplotlib import pyplot
import seaborn as sns

import os
import pandas as pd
import numpy as np

from collections import Counter
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import AllChem

class Visualization():
    def __init__(self,type_cluster,scaffold_type, generator_name, receptor_name, save_options):
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
        self.receptor_name = receptor_name
        self.save_options = save_options
        self.results = pd.DataFrame()



    def control_folder_existing(self):
        '''Save function to folder 'img'. If the folder doesn't exist -> the folder will be created'''
        print("Control folder existing")
        path_to_folder = 'img'
        

        #seconf to save barplots:
        if not os.path.exists(f"{path_to_folder}/plots/{self.receptor_name}/{self.type_cluster}/{self.scaffold_type}_scaffolds/{self.generator_name}/"):
            os.makedirs(f"{path_to_folder}/plots/{self.receptor_name}/{self.type_cluster}/{self.scaffold_type}_scaffolds/{self.generator_name}/")


        #third for save the most common
        if not os.path.exists(f"{path_to_folder}/the_most_common_scaffolds/{self.receptor_name}/{self.type_cluster}/{self.scaffold_type}_scaffolds/{self.generator_name}/"):
            os.makedirs(f"{path_to_folder}/the_most_common_scaffolds/{self.receptor_name}/{self.type_cluster}/{self.scaffold_type}_scaffolds/{self.generator_name}/")

    def save_dataframe_image(self, df, path_to_save, name, size_1, size_2):
        '''Save dataframe like image and then it's can be used in report'''
        shtml='''<html>
        <head>
        <style>
        .dataframe {
            background-color: white;
            color: black;
            padding: 200px;
            overflow-wrap: break-word;
            border-collapse: collapse;
        }
        </style>
        </head>
        <body>
        '''

        ehtml='''</body>
        </html>
        '''

        print("EXPORT_using_script")
        hti = Html2Image(output_path=path_to_save)

        ht=df.to_html(escape=False,index=False,border=1)
        ht=ht.replace("text-align: right","text-align: center")
        ht=shtml+ht+ehtml


        hti.screenshot(html_str=ht, save_as=name,size=(size_1, size_2))



class Visualization_mean():
    def __init__(self,type_cluster,scaffold_type, generators_name_list, receptor_name, save_options):
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
        self.generators_name_list = generators_name_list
        self.receptor_name = receptor_name
        self.save_options = save_options
        self.results = pd.DataFrame()



    def control_folder_existing(self):
        '''Save function to folder 'img'. If the folder doesn't exist -> the folder will be created'''
        print("Control folder existing")
        path_to_folder = 'img'
        

        #seconf to save barplots:
        if not os.path.exists(f"{path_to_folder}/plots/{self.receptor_name}/{self.type_cluster}/{self.scaffold_type}_scaffolds/compare/"):
            os.makedirs(f"{path_to_folder}/plots/{self.receptor_name}/{self.type_cluster}/{self.scaffold_type}_scaffolds/compare/")


    def save_dataframe_image(self, df, path_to_save, name, size_1, size_2):
        '''Save dataframe like image and then it's can be used in report'''
        shtml='''<html>
        <head>
        <style>
        .dataframe {
            background-color: white;
            color: black;
            padding: 200px;
            overflow-wrap: break-word;
            border-collapse: collapse;
        }
        </style>
        </head>
        <body>
        '''

        ehtml='''</body>
        </html>
        '''

        print("EXPORT_using_script")
        hti = Html2Image(output_path=path_to_save)

        ht=df.to_html(escape=False,index=False,border=1)
        ht=ht.replace("text-align: right","text-align: center")
        ht=shtml+ht+ehtml


        hti.screenshot(html_str=ht, save_as=name,size=(size_1, size_2))