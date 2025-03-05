import pandas as pd
import numpy as np
import time
from rdkit import Chem
from rdkit.Chem import AllChem
import umap
import matplotlib.pyplot as plt
import math
import os
from scipy.spatial.distance import jaccard
from scipy.spatial.distance import pdist, squareform



# Function to convert SMILES to Morgan Fingerprint
def smiles_to_morgan(smiles, radius=3, nbits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        return np.array(fingerprint)
    else:
        return np.zeros(nbits)



class Umap_class:
    def __init__(self, generators, type_cluster, number, receptor, labels_dict, num_cpus = 1):
        self.generators = generators
        self.type_cluster = type_cluster
        self.number = number
        self.receptor = receptor
        self.labels_dict = labels_dict
        self.num_cpus = num_cpus

        self.distance_matrix = None
        self.fingerprints = None


    def jaccard_distance(self, i, j):
        return jaccard(self.fingerprints[i], self.fingerprints[j])
    
    def process_for_plot_umap(self):
        """
        Process SMILES datasets, compute fingerprints, and generate UMAP visualization.
        """
        print(f"Processing dataset {self.number}...")
        before = time.perf_counter()

        # Load datasets
        datasets = []
        for i, generator in enumerate(self.generators):
            df = pd.read_csv(f"data/{generator}", header=None)
            data = df.sample(n=2500, random_state=42)[0].tolist() if i > 1 else df[0].tolist()
            datasets.append(data)

        # Prepare SMILES and labels
        combined_smiles = []
        labels = []
        for i, data in enumerate(datasets):
            for smiles in data:
                combined_smiles.append(smiles)
                labels.append(self.labels_dict[i])

        # Generate fingerprints
        print("Creating fingerprints...")
        self.fingerprints = np.array([smiles_to_morgan(smiles) for smiles in combined_smiles])
        print("Finished creating fingerprints")
        print('Len of fingeprints: ', len(self.fingerprints))

        # Compute Jaccard distance matrix using nested loops
        print("Calculating distance matrix...")
        num_samples = len(self.fingerprints)
        self.distance_matrix = np.zeros((num_samples, num_samples))
        
        chunks = 1000
        for i, fps1 in enumerate(self.fingerprints):
            if i % chunks == 0:
                    print("I: ", i)
            for j in range(i + 1, num_samples):
                value = self.jaccard_distance(i, j)
                self.distance_matrix[i][j] = value
                self.distance_matrix[j][i] = value

        # Apply UMAP
        umap_model = umap.UMAP(metric='precomputed', n_components=2, random_state=42)
        umap_results = umap_model.fit_transform(self.distance_matrix)

        # Store results
        umap_df = pd.DataFrame(umap_results, columns=['UMAP1', 'UMAP2'])
        umap_df['set_label'] = labels

        folder = f"data/results/{self.receptor}/UMAP"
        os.makedirs(folder, exist_ok=True)
        umap_df.to_csv(f"{folder}/umap_results_{self.type_cluster}_{self.number}.csv", index=False)

        after = time.perf_counter()
        print(f"Finished processing dataset {self.number} in {after - before:.2f} seconds.")




def plot_umap_results(type_cluster, number, receptor, unique_labels, name_save):
    """
    Function to visualize UMAP results for a given cluster type and dataset number.
    
    Arguments:
    type_cluster -- Type of clustering (default is 'sim').
    number -- Dataset number to process (default is 0).
    """
    # Load the UMAP results
    umap_df = pd.read_csv(f"data/results/{receptor}/UMAP/umap_results_{type_cluster}_{number}.csv")

    # Define color map
    colors = plt.get_cmap('tab10')

    # Create a plot with subplots
    fig, axes = plt.subplots( math.ceil(len(unique_labels)/2), 2, figsize=(30, math.ceil(len(unique_labels)/2)*7))

    # Loop over each axis and the corresponding labels
    for ax, labels in zip(axes.flatten(), unique_labels):
        alpha = 1
        for i, label in enumerate(labels):
            # Adjust color and alpha for 'IS' and 'RS' labels
            if label == 'IS':
                alpha = 1
                colors_ = 'gray'
            elif label == 'RS':
                alpha = 1
                colors_ = 'black'
            else:
                colors_ = colors(i)
            
            # Filter data for the current label
            subset = umap_df[umap_df['set_label'] == label]
            
            # Plot the scatter plot for each label
            ax.scatter(subset['UMAP1'], subset['UMAP2'], 
                       label=label, 
                       color=colors_, 
                       alpha=alpha) 

            # Adjust alpha based on the number of labels
            if len(labels) == 2:
                alpha -= 0.3
            else:
                alpha -= 0.15

        # Add title and axis labels
        title = ' vs. '.join(labels)
        ax.set_title(f'{title}', fontsize=18)  # Set larger title font size
        ax.set_xlabel('UMAP Component 1', fontsize=16)  # Set larger x-axis label font size
        ax.set_ylabel('UMAP Component 2', fontsize=16)  # Set larger y-axis label font size
        ax.legend(fontsize=12)
        ax.grid(True)
    plt.tight_layout()
    fig.suptitle(f'UMAP Visualizations for {type_cluster}_{number}', fontsize=20)
    fig.subplots_adjust(top=0.90)  # Můžeš experimentovat s hodnotou (např. 0.9, 0.95)
    
    # Save the plot as an image
    #plt.subplots_adjust(top=0.85) 

    folder = f"img/UMAP/{receptor}/double_plots/"
    os.makedirs(folder, exist_ok=True)

    plt.savefig(f'{folder}/{name_save}.png', format="png")
    # Display the plot
    plt.show()


def plot_umap_single_results(type_cluster, number, receptor, unique_labels, name_save):
    """
    Function to visualize UMAP single results for a given cluster type and dataset number.
    
    Arguments:
    type_cluster -- Type of clustering (default is 'sim').
    number -- Dataset number to process (default is 0).
    """
    # Load the UMAP results
    umap_df = pd.read_csv(f"data/results/{receptor}/UMAP/umap_results_{type_cluster}_{number}.csv")

    # Define color map
    colors = plt.get_cmap('tab10')

    # Create a plot with subplots
    fig, ax = plt.subplots( 1, 1, figsize=(12,8))

    # Loop over each axis and the corresponding labels
    
    # Loop over each axis and the corresponding labels
    for labels in unique_labels:
        alpha = 1
        for i, label in enumerate(labels):
            # Adjust color and alpha for 'IS' and 'RS' labels
            if label == 'IS':
                alpha = 1
                colors_ = 'gray'
            elif label == 'RS':
                alpha = 1
                colors_ = 'black'
            else:
                colors_ = colors(i)
            
            # Filter data for the current label
            subset = umap_df[umap_df['set_label'] == label]
            
            # Plot the scatter plot for each label
            ax.scatter(subset['UMAP1'], subset['UMAP2'], 
                       label=label, 
                       color=colors_, 
                       alpha=alpha) 

            # Adjust alpha based on the number of labels
            if len(labels) == 2:
                alpha -= 0.3
            else:
                alpha -= 0.15

        # Add title and axis labels
        title = ' vs. '.join(labels)
        ax.set_title(f'{title}', fontsize=18)  # Set larger title font size
        ax.set_xlabel('UMAP Component 1', fontsize=16)  # Set larger x-axis label font size
        ax.set_ylabel('UMAP Component 2', fontsize=16)  # Set larger y-axis label font size
        ax.legend(fontsize=12)
        ax.grid(True)
    
    fig.suptitle(f'UMAP Visualizations for {type_cluster}_{number} for {receptor}', fontsize=20)
    fig.subplots_adjust(top=0.90)  # Můžeš experimentovat s hodnotou (např. 0.9, 0.95)
    
    # Save the plot as an image
    #plt.subplots_adjust(top=0.85) 
    plt.tight_layout()
    folder = f"img/UMAP/{receptor}/single_plot/"
    os.makedirs(folder, exist_ok=True)
    
    plt.savefig(f'{folder}/{name_save}.png', format="png")
    # Display the plot
    plt.show()

