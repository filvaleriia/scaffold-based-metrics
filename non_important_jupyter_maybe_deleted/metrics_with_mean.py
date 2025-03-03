"""
This script calculates key metrics (TUPOR, SESY, ASER) for molecular scaffold analysis.
"""

import os
import pandas as pd
import numpy as np
from multiprocessing import Pool

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles, MakeScaffoldGeneric

import argparse



def add_columns_same_like_input_function(df_generated, test_set):
    """
    Create a DataFrame with three columns:
      - Column 0: Unique scaffolds in the Recall Set (test_set)
      - Column 1: Count of occurrences of each scaffold in df_generated
      - Column 2: Indicator (1 if count > 0, else 0)
    """
    df_generated = pd.DataFrame(df_generated)
    result = pd.DataFrame({'scaffold': test_set})
    counts = df_generated[0].value_counts()
    result['count'] = result['scaffold'].apply(lambda x: counts[x] if x in counts.index else 0)
    result['unique_indicator'] = result['count'].apply(lambda x: 1 if x > 0 else 0)
    return result



class Tc:
    """Class to compute Tanimoto similarity between two sets of molecules."""
    def __init__(self, df1, df2):
        self.df1_recall = df1
        self.df2_output = df2

    def find_the_most_close_molecule(self, fps_recall):
        """Find the most similar molecule in the output set."""
        return max(DataStructs.FingerprintSimilarity(fps_recall, fps_output) for fps_output in self.df2_output)



class Metrics:
    def __init__(self, type_cluster: str, type_scaffold: str, generator_name: str, receptor: str, 
                 save_options: bool, cal_median: bool = False, num_cpus: int = 1):
        self.type_cluster = type_cluster
        self.type_scaffold = type_scaffold
        self.generator_name = generator_name
        
        self.receptor = receptor
        self.save_options = save_options
        self.cal_meadian = cal_median
        self.num_cpus = num_cpus

        self.number_of_calculation = None
        self.output_set = None
        self.recall_set = None
        self.output_set_scaffolds = None
        self.recall_set_scaffolds = None
        self.unique_output_set = None
        self.unique_recall_set = None
        self.count_metrics = None
        self.results = pd.DataFrame()


    def convert_to_scaffold_(self, smiles_str):
        """
        Convert a SMILES string to its scaffold.
        For 'csk', return the generic Murcko scaffold.
        For 'murcko', return the Murcko scaffold.
        """
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is None:
                return None
            if self.type_scaffold == 'csk':
                scaffold = MurckoScaffoldSmiles(Chem.MolToSmiles(MakeScaffoldGeneric(mol)))
            elif self.type_scaffold == 'murcko':
                scaffold = MurckoScaffoldSmiles(Chem.MolToSmiles(mol))
            else:
                return None
            return scaffold if scaffold != '' else None
        except Exception:
            return None


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


    def scaffolds_unique(self, output_set, recall_set):
        """
        Return unique scaffolds from the output and recall sets.
        """
        return output_set[0].unique(), recall_set[0].unique()


    def main_function_return(self, output_set, recall_set):
        """
        Calculate all metrics.
        """
        # Convert compounds in recall set to scaffolds using multiprocessing.
        with Pool(processes=self.num_cpus) as pool:
            recall_results = pool.map(self.convert_to_scaffold_, recall_set[0])
        self.recall_set_scaffolds = pd.DataFrame(recall_results).dropna()

        # Convert compounds in output set to scaffolds.
        with Pool(processes=self.num_cpus) as pool:
            output_results = pool.map(self.convert_to_scaffold_, output_set[0])
        self.output_set_scaffolds = pd.DataFrame(output_results).dropna()

        # Get unique scaffolds for output and recall sets.
        self.unique_output_set, self.unique_recall_set = self.scaffolds_unique(self.output_set_scaffolds, self.recall_set_scaffolds)

        # Calculate occurrence of each scaffold in the output set.
        df = add_columns_same_like_input_function(self.output_set_scaffolds, self.unique_recall_set)
        self.count_metrics = df.copy()
        self.count_metrics.columns = ['unique_scaffold_recall', 'count_of_occurrence', 'unique_indicator']

        # Calculate metrics.
        USo = len(self.unique_output_set)
        SSo = len(self.output_set_scaffolds)
        CwASo = self.count_metrics['count_of_occurrence'].sum()
        try:
            UASo = self.count_metrics['unique_indicator'].value_counts()[1]
        except Exception:
            UASo = 0

        SSr = len(recall_set)
        UASr = len(self.count_metrics)

        try:
            tupor = UASo / UASr
            tupor_text = f"{UASo}/{UASr}"
        except Exception:
            tupor = 0
            tupor_text = f"0/{len(df)}"

        # Calculate additional metrics.
        other_scaffolds_unique = [scf for scf in self.unique_output_set if scf not in self.recall_set_scaffolds[0].tolist()]
        other_scaffolds = [scf for scf in self.output_set_scaffolds[0].tolist() if scf not in self.recall_set_scaffolds[0].tolist()]

        SESY = USo / SSo if SSo > 0 else 0
        ASER = CwASo / len(other_scaffolds) if len(other_scaffolds) > 0 else 0

        return self.type_cluster, USo, SSo, tupor_text, tupor, SESY, ASER, CwASo


    def save_function(self):
        """
        Save the calculated metrics to files under 'data/results'.
        """
        folder = f"data/results/{self.receptor}/{self.type_scaffold}_scaffolds/{self.type_cluster}/{self.generator_name}"
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        self.count_metrics.to_csv(f"{folder}/count_of_occurrence_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)
        self.results.to_csv(f"{folder}/metrics_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", index=False)
        self.output_set_scaffolds.to_csv(f"{folder}/scaffolds_of_output_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)
        self.recall_set_scaffolds.to_csv(f"{folder}/scaffolds_of_recall_set_cluster_{self.number_of_calculation}_{self.type_cluster}_{self.generator_name}.csv", header=False, index=False)


    def average_value(self, numbers):
        """
        For each number in 'numbers', load the corresponding CSV file,
        combine them into one DataFrame, compute the mean of numeric columns,
        append the mean row to the DataFrame, and save the results.
        """
        base_path = f"data/results/{self.receptor}/{self.type_scaffold}_scaffolds/{self.type_cluster}/{self.generator_name}/"

        # Build file paths for each cluster number.
        file_paths = {x: f"{base_path}metrics_cluster_{x}_{self.type_cluster}_{self.generator_name}.csv" for x in numbers}


        # Load all CSV files and combine them.
        combined_df = pd.concat([pd.read_csv(path) for path in file_paths.values()], ignore_index=True)

        # Calculate the mean for numeric columns.
        mean_values = combined_df.mean(numeric_only=True)

        # Create a mean row with the required format.
        # Use '-' as a placeholder for non-numeric columns.
        if self.cal_meadian:
            mean_row = [
                f"{self.generator_name}_mean",         # Name
                self.type_cluster,                # Type cluster
                self.type_scaffold,               # Scaffold type
                mean_values.get('USo', np.nan),
                mean_values.get('SSo', np.nan),
                '-',                         # Placeholder column
                mean_values.get('TUPOR', np.nan),
                mean_values.get('SESY', np.nan),
                mean_values.get('ASER', np.nan),
                mean_values.get('CwASo', np.nan),
                mean_values.get('median_OS_IS', np.nan),
            ]
        else:
            mean_row = [
                f"{self.generator_name}_mean",         # Name
                self.type_cluster,                # Type cluster
                self.type_scaffold,               # Scaffold type
                mean_values.get('USo', np.nan),
                mean_values.get('SSo', np.nan),
                '-',                         # Placeholder column
                mean_values.get('TUPOR', np.nan),
                mean_values.get('SESY', np.nan),
                mean_values.get('ASER', np.nan),
                mean_values.get('CwASo', np.nan)
            ]

        # Append the mean row to the combined DataFrame.
        combined_df.loc[len(combined_df)] = mean_row

        # Round numeric values to 7 decimals.
        combined_df = combined_df.round(7)

        # Create a formatted copy for specific columns.
        formatted_df = combined_df.copy()
        for col in ['SSo', 'USo', 'CwASo']:
            formatted_df[col] = formatted_df[col].apply(lambda x: "{:,}".format(x) if pd.notnull(x) else x)

        # Display the complete DataFrame and the mean row.

        mean_df = combined_df.tail(1)


        # Save the DataFrames.
        formatted_df.to_csv(f"{base_path}df_all_clusters_with_mean_with_coma.csv", index=False)
        combined_df.to_csv(f"{base_path}df_all_clusters_with_mean.csv", index=False)
        formatted_df.to_html(f"{base_path}df_all_clusters_with_mean.html", index=False)
        mean_df.to_csv(f"{base_path}{self.generator_name}_mean_{self.type_scaffold}_{self.type_cluster}.csv", index=False)

        return combined_df
    

    def KL_divergence_prep_mutliprocess(self, df1, df2):
        """
        Prepare data and compute similarity metrics (Tanimoto coefficient) between two sets using Morgan fingerprints.

        Args:
            df1 (pd.DataFrame): First DataFrame containing data. Each entry is a SMILES string representing a compound in the recall set.
            df2 (pd.DataFrame): Second DataFrame containing data. Each entry is a SMILES string representing a compound in the output set.

        Returns:
            list: List of highest Tanimoto coefficients between the two sets of compounds.
        """
        # Generate Morgan fingerprints for compounds in df1 and df2
        fps1 = []
        df1_copy = df1[0].tolist()  # Extracting SMILES from the first column of df1
        for x in df1_copy:
            try:
                fps1.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 3, nBits=2048))
            except:
                print("EXCEPT fps1", x)

        fps2 = []
        df2_copy = df2[0].tolist()  # Extracting SMILES from the first column of df2
        for x in df2_copy:
            try:
                fps2.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 3, nBits=2048))
            except:
                print("EXCEPT fps2", x)

        mt = Tc(fps1, fps2)

        # Use multiprocessing to compute the most similar molecules using Tanimoto coefficient
        with Pool(processes=self.num_cpus) as pool:
            tc_arr = pool.map(mt.find_the_most_close_molecule, fps1)

        return tc_arr


    def median_OS_IS(self, num):
        """
        Calculate the median of the Tanimoto coefficients between scaffolds of the output and input sets.

        Args:
            num (int): A number used in file naming, related to the cluster type or dataset index.

        Returns:
            float: The median of the Tanimoto coefficients for scaffold-based similarity between the output and input sets.
        """
        print("MEDIAN")
        input_set_path = f'data/input_recall_sets/{self.receptor}/cIS_{self.receptor}_{self.type_cluster}_{num}.csv'

        # Read the input file containing SMILES strings for the recall set
        with open(input_set_path, 'r') as f:
            input_set = f.read().splitlines()

        input_set_df = pd.DataFrame(input_set)

        # Use multiprocessing to convert SMILES to scaffolds
        with Pool(processes=self.num_cpus) as pool:
            input_results = pool.map(self.convert_to_scaffold_, input_set_df[0])

        input_set_scaffolds = pd.DataFrame(input_results).dropna()

        # Split the output set into chunks for parallel processing
        n_splits = 20
        len_df_output = len(self.output_set_scaffolds)
        print(f"Total length: {len_df_output}")

        chunk_size = len_df_output // n_splits
        chunks = [self.output_set_scaffolds.iloc[i * chunk_size: (i + 1) * chunk_size] for i in range(n_splits)]

        # Add remaining elements to the last chunk if necessary
        if len_df_output % n_splits != 0:
            chunks[-1] = pd.concat([chunks[-1], self.output_set_scaffolds.iloc[n_splits * chunk_size:]])

        # Calculate Tanimoto coefficients for each chunk sequentially
        tc_arr_all = []
        for i, chunk in enumerate(chunks):
            print(f"Processing chunk {i+1}/{n_splits}...")
            tc_arr = self.KL_divergence_prep_mutliprocess(chunk, input_set_scaffolds)
            tc_arr_all += tc_arr

        # Convert the results to a DataFrame and save them
        dff = pd.DataFrame(tc_arr_all)
        base_path = f"data/results/{self.receptor}/{self.type_scaffold}_scaffolds/{self.type_cluster}/{self.generator_name}/"
        os.makedirs(base_path, exist_ok=True)  # Create directories if they do not exist
        dff.to_csv(f"{base_path}high_TC_OS_IS_{self.type_scaffold}_{self.generator_name}_{self.type_cluster}_{num}.csv", header=None, index=False)

        # Return the median of the Tanimoto coefficients
        return np.median(dff)


    def calculate_metrics(self):
        """
        Calculate all metrics and return a DataFrame with the results.
        """
        numbers_of_calcs = []
        for number in range(5):
            print("NUMBER: ", number)
            self.number_of_calculation = number
            output_file_path = f"data/output_sets/{self.receptor}/{self.generator_name}/cOS_{self.generator_name}_{self.type_cluster}_{self.number_of_calculation}_one_column.csv"
            recall_file_path = f"data/input_recall_sets/{self.receptor}/cRS_{self.receptor}_{self.type_cluster}_{self.number_of_calculation}.csv"
            if os.path.exists(output_file_path):
                self.load(output_file_path, recall_file_path)
                res = self.main_function_return(self.output_set, self.recall_set)
                results = pd.DataFrame(columns=['type_cluster', 'USo', 'SSo', 'TUPOR_', 'TUPOR', 'SESY', 'ASER', 'CwASo'])
                results.loc[len(results)] = res
                results.insert(0, 'name', f"{self.generator_name}_{self.number_of_calculation}")
                results.insert(2, 'scaffold', self.type_scaffold)
                if self.cal_meadian:
                    median = self.median_OS_IS(self.number_of_calculation)
                    results.insert(len(results.columns), 'median_OS_IS', median)

                self.results = results

                if self.save_options:
                    self.save_function()
                numbers_of_calcs.append(self.number_of_calculation)
            else:
                print(f"Path for split {self.number_of_calculation} doesn't exists")
        
        
        result = self.average_value(numbers_of_calcs)
        
        
        return result
    

def main():
    parser = argparse.ArgumentParser(description='Compute and visualize recall metric.')
     # Required arguments
    parser.add_argument('--type_cluster', type=str, required=True, help='Type of clustering (dis/sim)')
    parser.add_argument('--type_scaffold', type=str, required=True, help='Type of scaffold')
    parser.add_argument('--generator', type=str, required=True, help='Generator name')
    parser.add_argument('--receptor', type=str, required=True, help='Receptor name')

    # Optional arguments with default values
    parser.add_argument('--save_option', type=bool, default=True, required=False, help='If set, results will be saved')
    parser.add_argument('--cal_median', type=bool, default=False, required=False, help='If set, median between OS and IS will be calculated')
    parser.add_argument('--ncpus', type=bool, default=1, required=False, help='Number of CPUs to use for parallel processing')

    
    args = parser.parse_args()
    print(args.cal_median)
    
    mt = Metrics(args.type_cluster, args.type_scaffold, args.generator, args.receptor, args.save_option, args.cal_median, args.ncpus)     
    result = mt.calculate_metrics()


if __name__ == "__main__":
    main()

