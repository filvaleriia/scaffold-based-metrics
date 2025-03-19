import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
from pathlib import Path
import argparse
import datetime
import sys
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate
from rdkit import DataStructs
from itertools import combinations
import ast




def tanimoto_similarity(vec1, vec2):
    # Vypočítá průnik a sjednocení
    intersection = np.count_nonzero(np.logical_and(vec1, vec2))  # Počet společných 1's
    union = np.count_nonzero(np.logical_or(vec1, vec2))  # Počet všech 1's v alespoň jednom vektoru
    
    # Vypočítá Tanimoto podobnost
    if union == 0:
        return 1.0  # Pokud jsou oba vektory prázdné (jen 0's), podobnost bude 1 (identické)
    else:
        return intersection / union

def create_matching_dataframe(recall_fps, output_fps):
    data = []
    data_fp = []
    # Převedeme recall_fps a output_fps na numpy array pro efektivní operace

    print(recall_fps)

    recall_fps = [np.array([int(bit) for bit in fp]) for fp in recall_fps[0].tolist()]
    output_fps = [np.array([int(bit) for bit in fp]) for fp in output_fps[0].tolist()]

    # Urychlení pomocí vektorových operací
    print("RECAL FP LEN: ", len(recall_fps))
    num = 0
    for recall_fp in recall_fps:
        print(recall_fp)


        # Počítání podobnosti Tanimoto pro všechny output_fps
        match_count = np.sum([tanimoto_similarity(recall_fp, output_fp) == 1.0 for output_fp in output_fps])
        
        # Pokud je počet shod > 0, přidáme řádek do seznamu
        data.append({
            'label': f'FP-{num}',
            'UAFo': 1 if match_count > 0 else 0,
            'CwAFo': match_count
        })
        data_fp.append({
            'label': f'FP-{num}',
            'recall_fingerprint': recall_fp,

        })
        num += 1
    print(data)

    return data


def average_tanimoto_diversity(fingerprints):
    print('START CALCULATE DIVERSITY')
    num_pairs = 0
    similarity_sum = 0
    for fp1, fp2 in combinations(fingerprints, 2):
        similarity_sum += DataStructs.FingerprintSimilarity(fp1, fp2)
        num_pairs += 1
    return similarity_sum / num_pairs if num_pairs > 0 else 0.0


class Metrics_phfp:
    def __init__(self, type_cluster: str, type_phfp: str, generator_name: str, receptor: str):
        self.type_cluster = type_cluster
        self.type_phfp = type_phfp
        self.generator_name = generator_name
        
        self.receptor = receptor
        self.number_of_calculation = None
        self.output_set = None
        self.recall_set = None
        self.output_set_phfp = None
        self.recall_set_phfp = None
        self.unique_output_set = None
        self.unique_recall_set = None
        self.count_metrics = None
        self.results = pd.DataFrame()


    def load(self, filepath_output_set, filepath_recall_set):
        self.output_set_phfp = pd.read_csv(filepath_output_set, header = None)[:10]
        self.recall_set_phfp = pd.read_csv(filepath_recall_set, header = None)[:10]
        print('ORIGINAL RECALL LEN: ', len(self.recall_set_phfp))
        self.recall_set_phfp = self.recall_set_phfp.drop_duplicates(keep='first').reset_index(drop=True)
        print('UNIQUE RECALL LEN: ', len(self.recall_set_phfp))


        
    def calculate_metrics(self):
        df = create_matching_dataframe(self.recall_set_phfp, self.output_set_phfp)
        self.count_metrics = pd.DataFrame(df)
        print("END CREATE MATCHING")
        # Calculate metrics.
        print(self.count_metrics['CwAFo'])
        UFo =  len(list({tuple(arr): arr for arr in self.output_set_phfp[0].tolist()}.values()))
        SSo = len(self.output_set_phfp)
        CwAFo = self.count_metrics['CwAFo'].sum()
        
        UAFo = self.count_metrics['UAFo'].sum()

        UAFr = len(self.count_metrics)

       
        tupor = UAFo / UAFr
        tupor_text = f"{UAFo}/{UAFr}"

        SESY = UFo / SSo 
        ASER = CwAFo / (len(self.output_set_phfp) - CwAFo)
        ACR = UAFo/(UFo-UAFo)
        print("END CALCULATE METRICS")

        #diversity = average_tanimoto_diversity(self.output_set_phfp)
        diversity = 0
        print('DIVERSITY END')

        results = pd.DataFrame(columns=['type_cluster', 'UFo', 'SSo', 'TUPOR_', 'TUPOR', 'SESY', 'ASER', 'CwAFo', 'ACR', 'DIV'])
        results.loc[len(results)] = [self.type_cluster, UFo, SSo, tupor_text, tupor, SESY, ASER, CwAFo, ACR, diversity]
        results.insert(0, 'name', f"{self.generator_name}_{self.number_of_calculation}")
        results.insert(2, 'phfp', self.type_phfp)

        main_dir = Path(__file__).resolve().parents[1]
        folder = f"{main_dir}/data/results_phram_fp/{self.receptor}/{self.type_phfp}/{self.type_cluster}/{self.generator_name}/"
        self.count_metrics.to_csv(f"{folder}/count_of_occurrence_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)
        results.to_csv(f"{folder}/metrics_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)

        print(results)


    def average_value(self, numbers):
        """
        For each number in 'numbers', load the corresponding CSV file,
        combine them into one DataFrame, compute the mean of numeric columns,
        append the mean row to the DataFrame, and save the results.
        """
        main_dir = Path(__file__).resolve().parents[1] 
        base_path = f"{main_dir}/data/results_phram_fp/{self.receptor}/{self.type_phfp}/{self.type_cluster}/{self.generator_name}/"

        # Build file paths for each cluster number.
        file_paths = {x: f"{base_path}metrics_cluster_{x}_{self.type_cluster}_{self.generator_name}.csv" for x in numbers}


        # Load all CSV files and combine them.
        combined_df = pd.concat([pd.read_csv(path) for path in file_paths.values()], ignore_index=True)

        # Calculate the mean for numeric columns.
        mean_values = combined_df.mean(numeric_only=True)

        # Create a mean row with the required format.
        
        mean_row = [
                f"{self.generator_name}_mean",         # Name
                self.type_cluster,                # Type cluster
                self.type_phfp,               # Scaffold type
                mean_values.get('UFo', np.nan),
                mean_values.get('SSo', np.nan),
                '-',                         # Placeholder column
                mean_values.get('TUPOR', np.nan),
                mean_values.get('SESY', np.nan),
                mean_values.get('ASER', np.nan),
                mean_values.get('CwAFo', np.nan),
                mean_values.get('ACR', np.nan),
                mean_values.get('DIV', np.nan),
            ]

        # Append the mean row to the combined DataFrame.
        combined_df.loc[len(combined_df)] = mean_row

        # Round numeric values to 7 decimals.
        combined_df = combined_df.round(7)

        # Create a formatted copy for specific columns.
        formatted_df = combined_df.copy()
        for col in ['SSo', 'UFo', 'CwAFo']:
            formatted_df[col] = formatted_df[col].apply(lambda x: "{:,}".format(x) if pd.notnull(x) else x)

        # Display the complete DataFrame and the mean row.

        mean_df = combined_df.tail(1)


        # Save the DataFrames.
        formatted_df.to_csv(f"{base_path}df_all_clusters_with_mean_with_coma.csv", index=False)
        combined_df.to_csv(f"{base_path}df_all_clusters_with_mean.csv", index=False)
        mean_df.to_csv(f"{base_path}{self.generator_name}_mean_{self.type_phfp}_{self.type_cluster}.csv", index=False)

        return combined_df


    def calculate(self):
        main_dir = Path(__file__).resolve().parents[1] 
        base_path = f"{main_dir}/data/results_phram_fp/{self.receptor}/{self.type_phfp}/{self.type_cluster}/{self.generator_name}/"

        numbers = []
        for number in range(5):
            self.number_of_calculation = number
            output_file_path = f"{base_path}phfp_of_output_set_cluster_{number}_{self.type_cluster}_{self.generator_name}.csv"
            recall_file_path = f"{base_path}phfp_of_recall_set_cluster_{number}_{self.type_cluster}_{self.generator_name}.csv"
            print(output_file_path)
            if os.path.exists(output_file_path):
                print("EXIST")
                self.load(output_file_path, recall_file_path)
                numbers.append(number)
                self.calculate_metrics()
        if numbers:
            self.average_value(numbers)
        else:
            print('No data for calculation')



def main():
    parser = argparse.ArgumentParser(description='Compute and visualize recall metric.')
    parser.add_argument('--type_cluster', type=str, required=True, help='Type of clustering (dis/sim)')
    parser.add_argument('--type_phfp', type=str, required=True, help='Type of scaffold')
    parser.add_argument('--generator', type=str, required=True, help='Generator name')
    parser.add_argument('--receptor', type=str, required=True, help='Receptor name')

    args = parser.parse_args()
    mt = Metrics_phfp(args.type_cluster, args.type_phfp, args.generator, args.receptor)
    mt.calculate()


if __name__ == "__main__":
    main()
