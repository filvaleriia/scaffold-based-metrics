"""
This script calculates key metrics (TUPOR, SESY, ASER) for molecular scaffold analysis.
"""

import os
import subprocess
import pandas as pd
import numpy as np
from multiprocessing import Pool
from itertools import combinations

import rdkit
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate
from rdkit import DataStructs

import argparse
from pathlib import Path

import CDPL
import CDPL.Chem as CDPLChem
import CDPL.Pharm as CDPLPharm
import CDPL.Util as Util
import CDPL.Descr as Descr


def create_matching_dataframe(recall_fps, output_fps):
    data = []
    
    for recall_fp in recall_fps:
        # Počítání kolikrát je recall_fp nalezen v output_fps
        match_count = sum(1 for output_fp in output_fps if DataStructs.TanimotoSimilarity(recall_fp, output_fp) == 1.0)
        
        # Pokud je počet shod > 0, přidáme řádek do seznamu
        data.append({
            'recall_fingerprint': recall_fp.ToBitString(),
            'UAFo': 1 if match_count > 0 else 0,
            'CwAFo': match_count
        })
    
    # Vytvoření DataFrame
    df = pd.DataFrame(data)
    return df


def prepradet_for_2Df_rdkit():
    main_dir = Path(__file__).resolve().parents[1] 
    fdefName = f'{main_dir}/BaseFeatures.fdef'
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    
    sigFactory = SigFactory(featFactory,minPointCount=2,maxPointCount=3)
    sigFactory.SetBins([(0,2),(2,4),(4,6),(6,8)])
    sigFactory.Init()
    return sigFactory


def average_tanimoto_diversity(fingerprints):
    num_pairs = 0
    similarity_sum = 0
    for fp1, fp2 in combinations(fingerprints, 2):
        similarity_sum += DataStructs.FingerprintSimilarity(fp1, fp2)
        num_pairs += 1
    return similarity_sum / num_pairs if num_pairs > 0 else 0.0


class Metrics_phfp:
    def __init__(self, type_cluster: str, type_phfp: str, generator_name: str, receptor: str, number: int, num_cpus: int = 1):
        self.type_cluster = type_cluster
        self.type_phfp = type_phfp
        self.generator_name = generator_name
        
        self.receptor = receptor
        self.num_cpus = num_cpus
        self.number = number

        self.number_of_calculation = None
        self.output_set = None
        self.recall_set = None
        self.output_set_phfp = None
        self.recall_set_phfp = None
        self.unique_output_set = None
        self.unique_recall_set = None
        self.count_metrics = None
        self.results = pd.DataFrame()
        

    def convert_to_phfp_rdkit(self, mol):
        #print("MOL: ", mol)
        try:
            sigFactory = prepradet_for_2Df_rdkit()

            fp = Generate.Gen2DFingerprint(Chem.MolFromSmiles(mol), sigFactory)
            return np.array(fp)
        except Exception as e:
            pass

    def convert_to_phfp_cdpkit(self, mol: CDPLChem.Molecule, num_bits: int, bin_size: float) -> Util.BitSet:
        
        CDPLPharm.prepareForPharmacophoreGeneration(mol)  # prepare molecule for pharmacophore generation
        fp_gen = Descr.NPoint2DPharmacophoreFingerprintGenerator() # create 2D pham. fingerprint generator instance
        fp_gen.setBinSize(bin_size) # set feature distance bin size     
        fp = Util.BitSet()       # create fingerprint bitset
        fp.resize(num_bits)      # set desired fingerprint size
        fp_gen.generate(mol, fp) # generate the fingerprint

        return fp



    def load(self, filepath_output_set, filepath_recall_set):
        """
        Load the output and recall sets from the given file paths.
        """
        
        with open(filepath_output_set, 'r') as f:
            output_set = f.read().splitlines()
        self.output_set = pd.DataFrame(output_set)

        with open(filepath_recall_set, 'r') as f:
            recall_set = f.read().splitlines()
        self.recall_set = pd.DataFrame(recall_set)


    def calculate_phfp(self, output_set, recall_set):
        """
        Calculate all metrics.
        """
        print("MAIN FUNCTION")
        print("RECALL SET LEBN: ", len(recall_set))
        # Convert compounds in recall set to scaffolds using multiprocessing.

        if self.type_phfp == 'rdkit':
            with Pool(processes=self.num_cpus) as pool:
                recall_results = pool.map(self.convert_to_phfp_rdkit, recall_set[0].tolist())
                pool.close()
                pool.join()
            self.recall_set_phfp = [fp for fp in recall_results if fp is not None]

            # Get unique scaffolds for recall sets.
            self.unique_recall_set = list({tuple(arr): arr for arr in self.recall_set_phfp}.values())

            print('UNIQUE recall: ', len(self.unique_recall_set))
            print("RECALL CONVERT DONE")
            print('OUTPUT_SET_LEN: ', len(output_set))
            # Convert compounds in output set to scaffolds.

            res_all = []
            output_parts = np.array_split(output_set[0], 10)

            for part in output_parts:
                with Pool(processes=self.num_cpus) as pool:
                    # Zpracování aktuální části
                    print("LEN part: ", len(part))
                    first_results = pool.map(self.convert_to_phfp_rdkit, part)
                    pool.close()
                    pool.join()
                print('Done')
                # Přidání výsledků do hlavního seznamu
                res_all.extend(first_results)
                print("RES ALL LEN: ", len(res_all))

            


        self.output_set_phfp = [fp for fp in res_all if fp is not None]

        print("OUTPUT CONVERT DONE")

        #save output phfp
        main_dir = Path(__file__).resolve().parents[1]
        folder = f"{main_dir}/data/results_phram_fp/{self.receptor}/{self.type_phfp}/{self.type_cluster}/{self.generator_name}"
        if not os.path.exists(folder):
            os.makedirs(folder)
        output_set_phfp_df = pd.DataFrame(data = self.output_set_phfp)
        output_set_phfp_df.to_csv(f"{folder}/phfp_of_output_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)

        recall_set_phfp_df = pd.DataFrame(data = self.recall_set_phfp)
        recall_set_phfp_df.to_csv(f"{folder}/phfp_of_recall_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)




    def calculate_cdpkit(self, output_set_path, recall_set_path):
        """
        Calculate all metrics.
        """
        print("MAIN FUNCTION")

        # Convert compounds in recall set to scaffolds using multiprocessing.
        main_dir = Path(__file__).resolve().parents[1]
        folder = f"{main_dir}/data/results_phram_fp/{self.receptor}/{self.type_phfp}/{self.type_cluster}/{self.generator_name}"
        if not os.path.exists(folder):
            os.makedirs(folder)
        

        smi_recall_file = recall_set_path.replace('.csv', '.smi')
        subprocess.run(['cp', recall_set_path, smi_recall_file])
        #RECALLL
        reader = CDPLChem.MoleculeReader(smi_recall_file)
        mol = CDPLChem.BasicMolecule()
        out_file = open(f'{folder}/phfp_of_recall_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv', 'w')
        while reader.read(mol):
            try:
                fp = self.convert_to_phfp_cdpkit(mol, 4096, 1)
                out_file.write(str(fp))
                out_file.write('\n')
            except:
                pass

        #OUTPUT
        input_file = output_set_path
        temp_file = output_set_path.replace('column.', 'column_temp.')

        
        with open(input_file, 'r') as infile, open(temp_file, 'w') as outfile:
            for line in infile:
                # Nahradí 'si' správným 'Si' pouze jako prvek (ne součást slova)
                corrected_line = line.replace('si', 'Si')
                outfile.write(corrected_line)

        # Nahradíme původní soubor opraveným
        os.replace(temp_file, input_file)

        smi_output_file = output_set_path.replace('.csv', '.smi')
        subprocess.run(['cp', output_set_path, smi_output_file])
        reader = CDPLChem.MoleculeReader(smi_output_file)
        mol = CDPLChem.BasicMolecule()
        out_file = open(f'{folder}/phfp_of_output_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv', 'w')

        try: 
            while reader.read(mol):
                fp = self.convert_to_phfp_cdpkit(mol, 4096, 1)
                out_file.write(str(fp))
                out_file.write('\n')
        except (ValueError, TypeError) as e:
            print(f'Nepovedlo nacist vstup: {e}')
            print(f'pro output: {smi_output_file}')




    def main_function_return(self, output_set, recall_set):
        """
        Calculate all metrics.
        """
        print("MAIN FUNCTION")
        print("RECALL SET LEBN: ", len(recall_set))
        # Convert compounds in recall set to scaffolds using multiprocessing.
        with Pool(processes=self.num_cpus) as pool:
            recall_results = pool.map(self.convert_to_phfp_rdkit, recall_set[0].tolist())
        self.recall_set_phfp = [fp for fp in recall_results if fp is not None]

        # Get unique scaffolds for recall sets.
        self.unique_recall_set = list({tuple(arr): arr for arr in self.recall_set_phfp}.values())

        print('UNIQUE recall: ', len(self.unique_recall_set))
        print("RECALL CONVERT DONE")
        print('OUTPUT_SET_LEN: ', len(output_set))
        # Convert compounds in output set to scaffolds.

        res_all = []
        output_parts = np.array_split(output_set[0], 10)

        for part in output_parts:
            with Pool(processes=self.num_cpus) as pool:
                # Zpracování aktuální části
                print("LEN part: ", len(part))
                first_results = pool.map(self.convert_to_phfp_rdkit, part)
                pool.close()
                pool.join()
            

            print('Done')
            # Přidání výsledků do hlavního seznamu
            res_all.extend(first_results)
            print("RES ALL LEN: ", len(res_all))

        #count = 0
        #for mol in output_set[0].tolist():
        #    res = self.convert_to_phfp_rdkit(mol)
        #    res_all.append(res)
        #    count += 1
        #    if count % 1000 == 0:
        #        print(count)

        

        # Spojení všech výsledků a filtrování None hodnot
        #print(res_all)
        self.output_set_phfp = [fp for fp in res_all if fp is not None]

        print("OUTPUT CONVERT DONE")

        #save output phfp
        main_dir = Path(__file__).resolve().parents[1]
        folder = f"{main_dir}/data/results_phram_fp/{self.receptor}/{self.type_phfp}/{self.type_cluster}/{self.generator_name}"
        if not os.path.exists(folder):
            os.makedirs(folder)
        output_set_phfp_df = pd.DataFrame(data = self.output_set_phfp)
        output_set_phfp_df.to_csv(f"{folder}/phfp_of_output_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)

        recall_set_phfp_df = pd.DataFrame(data = self.recall_set_phfp)
        recall_set_phfp_df.to_csv(f"{folder}/phfp_of_recall_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)
        
        print('create matching')
        # Calculate occurrence of each scaffold in the output set.
        df = create_matching_dataframe(self.unique_recall_set, self.output_set_phfp)
        self.count_metrics = df.copy()

        print('end create matching')
        # Calculate metrics.
        UFo =  len(list({tuple(arr): arr for arr in self.output_set_phfp}.values()))

        SSo = len(self.output_set_phfp)

        CwAFo = self.count_metrics['CwAFo'].sum()

        UAFo = self.count_metrics['UAFo'].sum()
        
        UAFr = len(self.count_metrics)

        try:
            tupor_pharm = UAFo / UAFr
            tupor_pharm_text = f"{UAFo}/{UAFr}"
        except Exception:
            tupor_pharm = 0
            tupor_pharm_text = f"0/{len(df)}"

        SESY_pharm = UFo / SSo
        ASER_pharm = CwAFo / (SSo - CwAFo)
        ACR_pharm = UAFo / (UFo - UAFo)

        print('create diversity')
        diversity = average_tanimoto_diversity(self.output_set_phfp)

        return self.type_cluster, UFo, SSo, tupor_pharm_text, tupor_pharm, SESY_pharm, ASER_pharm, ACR_pharm, CwAFo, diversity



    def save_function(self):
        """
        Save the calculated metrics to files under 'data/results_pharm_fp'.
        """
        main_dir = Path(__file__).resolve().parents[1]
        folder = f"{main_dir}/data/results_phram_fp/{self.receptor}/{self.type_phfp}/{self.type_cluster}/{self.generator_name}"
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        self.count_metrics.to_csv(f"{folder}/count_of_occurrence_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)
        self.results.to_csv(f"{folder}/metrics_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)

        output_set_phfp_df = pd.DataFrame(data = self.output_set_phfp)
        output_set_phfp_df.to_csv(f"{folder}/phfp_of_output_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)

        recall_set_phfp_df = pd.DataFrame(data = self.recall_set_phfp)
        recall_set_phfp_df.to_csv(f"{folder}/phfp_of_recall_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)


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
        # Use '-' as a placeholder for non-numeric columns.
        
        mean_row = [
                f"{self.generator_name}_mean",         # Name
                self.type_cluster,                # Type cluster
                self.type_phfp,               # Scaffold type
                mean_values.get('UFo', np.nan),
                mean_values.get('SSo', np.nan),
                '-',                         # Placeholder column
                mean_values.get('TUPOR_pharm', np.nan),
                mean_values.get('SESY_pharm', np.nan),
                mean_values.get('ASER_pharm', np.nan),
                mean_values.get('ACR_pharm', np.nan),
                mean_values.get('CwAFo', np.nan),
                mean_values.get('diversita_pharm', np.nan)
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
    




    def calculate_metrics(self):
        """
        Calculate all metrics and return a DataFrame with the results.
        """
        main_dir = Path(__file__).resolve().parents[1] 
        numbers_of_calcs = []
        
        
        print("NUMBER: ", self.number)
        self.number_of_calculation = self.number
        output_file_path = f"{main_dir}/data/output_sets/{self.receptor}/{self.generator_name}/cOS_{self.generator_name}_{self.type_cluster}_{self.number_of_calculation}_one_column.csv"
        recall_file_path = f"{main_dir}/data/input_recall_sets/{self.receptor}/cRS_{self.receptor}_{self.type_cluster}_{self.number_of_calculation}.csv"
        print(output_file_path)
        if os.path.exists(output_file_path):
            print('EXIST')
            if self.type_phfp == 'rdkit':
                self.load(output_file_path, recall_file_path)
                self.calculate_phfp(self.output_set, self.recall_set)
            elif self.type_phfp == 'cdpkit':
                self.calculate_cdpkit(output_file_path, recall_file_path)
            
            '''
            results = pd.DataFrame(columns=['type_cluster', 'UFo', 'SSo', 'TUPOR_', 'TUPOR_pharm', 'SESY_pharm', 'ASER_pharm', 'ACR_pharm','CwAFo', 'diversita_pharm'])
            results.loc[len(results)] = res
            results.insert(0, 'name', f"{self.generator_name}_{self.number_of_calculation}")
            results.insert(2, 'phfp', self.type_phfp)
            self.results = results
            self.save_function()
            numbers_of_calcs.append(self.number_of_calculation)
            print(self.results)
        else:
            print(f"Path for cluster {self.number_of_calculation} doesn't exists")
    
    
        result = self.average_value(numbers_of_calcs)
        
        
        return result
        '''
    

def main():
    parser = argparse.ArgumentParser(description='Compute and visualize recall metric.')
     # Required arguments
    parser.add_argument('--type_cluster', type=str, required=True, help='Type of clustering (dis/sim)')
    parser.add_argument('--type_phfp', type=str, required=True, help='Type of scaffold')
    parser.add_argument('--generator', type=str, required=True, help='Generator name')
    parser.add_argument('--receptor', type=str, required=True, help='Receptor name')

    # Optional arguments with default values
    parser.add_argument('--ncpus', type=bool, default=1, required=False, help='Number of CPUs to use for parallel processing')

    
    args = parser.parse_args()
    print(args)
    
    mt = Metrics_phfp(args.type_cluster, args.type_phfp, args.generator, args.receptor, args.ncpus)     
    result = mt.calculate_metrics()


if __name__ == "__main__":
    main()

