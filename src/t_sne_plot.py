import pandas as pd
import numpy as np
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import math


# Function to convert SMILES to Morgan Fingerprint
def smiles_to_morgan(smiles, radius=3, nbits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        return np.array(fingerprint)
    else:
        return np.zeros(nbits)



def process_for_plot_tsne(generators, type_cluster, number, receptor, labels_dict):
    """
    Function to process SMILES datasets and plot t-SNE visualization for the specified generators and cluster type.

    Arguments:
    generators -- List of generator names (strings).
    type_cluster -- Type of clustering to use ('sim' or 'dis').
    numbers -- List of dataset numbers to process (default is [0, 1, 2, 3, 4]).
    """


    # Iterate over the dataset numbers
    
    print(f"Processing dataset {number}...")
    before = time.perf_counter()
    # List to hold datasets for all generators
    datasets = []
    # Load datasets for each generator and add them to the datasets list
    for i, generator in enumerate(generators):
        if i >1:
            datasets.append(pd.read_csv(f"data/{generator}", header=None) .sample(n=10000, random_state=42)[0].tolist())
        else:
            datasets.append(pd.read_csv(f"data/{generator}", header=None)[0].tolist())
    # Combine all SMILES data into one list
    combined_smiles = []
    labels = []
    for i, data in enumerate(datasets):
        combined_smiles.extend(data)
        labels.extend([labels_dict[i]] * len(data))  # Add the corresponding label for each dataset
    
    # Generate Morgan fingerprints for all SMILES
    fingerprints = np.array([smiles_to_morgan(smiles) for smiles in combined_smiles])
    
    # Apply t-SNE to reduce the dimensionality
    tsne = TSNE(n_components=2, random_state=42, perplexity=10)
    tsne_results = tsne.fit_transform(fingerprints)
    
    # Create a DataFrame to store t-SNE results
    tsne_df = pd.DataFrame(tsne_results, columns=['TSNE1', 'TSNE2'])
    tsne_df['set_label'] = labels  # Add the dataset labels
    
    # Save the t-SNE results to a CSV file
    tsne_df.to_csv(f"data/results/{receptor}/TSNE/tsne_results_{type_cluster}_{number}.csv", index=False)
    after = time.perf_counter()
    print(f"Finished processing dataset {number} in {after - before:.2f} seconds.")



def plot_tsne_results(type_cluster, number, receptor, unique_labels, name_save):
    """
    Function to visualize t-SNE results for a given cluster type and dataset number.
    
    Arguments:
    type_cluster -- Type of clustering (default is 'sim').
    number -- Dataset number to process (default is 0).
    """
    # Load the t-SNE results
    tsne_df = pd.read_csv(f"data/results/{receptor}/TSNE/tsne_results_{type_cluster}_{number}.csv")

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
            subset = tsne_df[tsne_df['set_label'] == label]
            
            # Plot the scatter plot for each label
            ax.scatter(subset['TSNE1'], subset['TSNE2'], 
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
        ax.set_xlabel('t-SNE Component 1', fontsize=16)  # Set larger x-axis label font size
        ax.set_ylabel('t-SNE Component 2', fontsize=16)  # Set larger y-axis label font size
        ax.legend(fontsize=12)
        ax.grid(True)
    plt.tight_layout()
    fig.suptitle(f't-SNE Visualizations for {type_cluster}_{number}', fontsize=20)
    fig.subplots_adjust(top=0.90)  # Můžeš experimentovat s hodnotou (např. 0.9, 0.95)
    
    # Save the plot as an image
    #plt.subplots_adjust(top=0.85) 
    plt.savefig(f'img/t-SNE/{receptor}/{name_save}.png', format="png")
    # Display the plot
    plt.show()


