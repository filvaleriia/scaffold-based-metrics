import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import MinMaxScaler


def preprocesing(type_cluster, type_scaffold, generators_name_list, receptor, ph4 = False):
    # Define path to data
    if ph4 == False:
        link = f"data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}"
    else:
        link = f"data/results_phram_fp/{receptor}/{type_scaffold}/{type_cluster}"
    link_mean = [f"{link}/{generator}/{generator}_mean_{type_scaffold}_{type_cluster}.csv" for generator in generators_name_list]
    
    # Load data
    df_list = [pd.read_csv(f) for f in link_mean]
    df = pd.concat(df_list, axis=0, ignore_index=True)

    scaler = MinMaxScaler()
    numeric_columns = df.select_dtypes(include=['number']).columns  # Select only numeric columns
    df[numeric_columns] = scaler.fit_transform(df[numeric_columns])  # Apply normalization
    df["name"] = df["name"].str.replace("_mean", "", regex=False)

    
    return df


def plot_heatmap(data, title='', name_save='',receptor = '', cmap='viridis', annotate=True, ph4 = False):
    ''' 
    Plots a single heatmap for the given data split.
    
    Args:
    - data (pd.DataFrame): The input data containing the columns to be visualized.
    - title (str): The title of the heatmap (default is empty).
    - cmap (str): The color map to be used for the heatmap (default is 'viridis').
    - annotate (bool): Whether to annotate the cells with their values (default is True).
    '''

    # Extract relevant columns (TUPOR, SESY, ASER, ASR) for visualization
    df = data[['TUPOR', 'SESY', 'ASER']]
    # Set the index of the dataframe to the 'name' attribute of the data
    df.index = data.name.tolist()

    y_labels = {

    "Molpher": "Molpher",
    "REINVENT": "REINVENT",
    "DrugEx_GT_epsilon_0.6" : "DrugEx_GT",
    "DrugEx_RNN_epsilon_0.6": "DrugEx_RNN",
    "GB_GA_mut_r_0.5": "GB_GA",
    "addcarbon": "AddCarbon"
    }
    df.rename(index=y_labels, inplace=True)

    # Create a figure for the heatmap with a specific size
    plt.figure(figsize=(10, 6))
    
    # Plot the heatmap using seaborn with optional annotations and custom color map
    sns.heatmap(df, annot=annotate, cmap=cmap,   annot_kws={"size": 17})
    
    # Set the title for the heatmap
    plt.title(title, fontsize=17 , pad=10)

    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    plt.tight_layout()
    # Save the plot as an SVG file
    if ph4:
        plt.savefig(f'img_pharm/heat_mapa/{receptor}/{name_save}.svg', format="svg")
        plt.savefig(f'img_pharm/heat_mapa/{receptor}/{name_save}.png', format="png")
    else:
        plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
        plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.png', format="png")
    # Display the heatmap
    plt.show()



def plot_all_subsets(subset_dict, title='', receptor = '', name_save = '', cmap='viridis', annotate=True, poradi = ''):
    '''
    Plots heatmaps for multiple subsets in a single figure.
    
    Args:
    - subset_dict (dict): Dictionary where keys are subset names (e.g., '', '_10k') 
                           and values are corresponding DataFrames to be visualized.
    - title (str): Title for the entire figure (default is empty).
    - cmap (str): Color map for the heatmaps (default is 'viridis').
    - annotate (bool): Whether to annotate the cells with their values (default is True).
    '''
    # Get the number of subsets in the dictionary
    num_subsets = len(subset_dict)
    
    # Create a subplot for each subset (1 row, num_subsets columns)
    fig, axes = plt.subplots(1, num_subsets, figsize=(num_subsets * 12, 12))
    
    # If there is only one subset, axes will not be iterable, so convert it to a list
    if num_subsets == 1:
        axes = [axes]
    
    # Iterate over each axis (subplot) and the corresponding subset data
    for ax, (subset_name, data) in zip(axes, subset_dict.items()):
        # Extract the relevant columns for the heatmap (TUPOR, SESY, ASER, ACR)
        df = data[['TUPOR', 'SESY', 'ASER']]
        
        # Set the index of the dataframe to the 'name' attribute of the data
        df.index = data.name.tolist()  # Using names as index
        
        # Plot the heatmap for the current subset
        sns.heatmap(df, annot=annotate, cmap=cmap, ax=ax,  annot_kws={"size": 30})

        
        # If the subset name is empty, label it as 'base'
        if subset_name == '':
            subset_name = 'Full OS'
        elif subset_name == '_62.5k':
            subset_name = '62,500'
        else:
            subset_name = subset_name.replace('_', '').replace('k', ',000')
        
        # Modify the y-axis labels for better readability by replacing certain substrings
        new_labels = [label.get_text().replace('_epsilon', '\n epsilon').replace('_mut_r', '\n mut_r').replace('addcarbon', 'AddCarbon') for label in ax.get_yticklabels()]
        new_labels = [label.replace('_62.5k', '').replace('_125k', '').replace('_250k', '').replace('_500k', '') for label in new_labels]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=30)
        ax.set_xticklabels(ax.get_xticklabels(), ha="center", fontsize=30)

        
        # Set the title for the current subplot to indicate the subset name
        if subset_name == 'Full OS':
            ax.set_title(f"{subset_name}",  fontsize=35, wrap=True)
        else:
            ax.set_title(f"{subset_name} subset",  fontsize=35, wrap=True)

    fig.text(
    0.005, 0.97, poradi,
    ha='left', va='top',
    fontsize=40
    )
    
    # Set the overall title for the figure
    fig.suptitle(f'{title}', fontsize=40)
    
    # Adjust layout to ensure titles and labels are well placed
    plt.tight_layout()
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.png', format="png")
    # Display the heatmap figure
    plt.show()



def plot_heatmap_base(subset_dict, subset_dict_data, title='', receptor = '', name_save = '', cmap='viridis', annotate=True, ph4 = False):
    '''
    Plots heatmaps for different subsets in a 2x2 grid, with each subset visualized in a separate subplot.
    
    Args:
    - subset_dict (dict): A dictionary containing subset names as keys (e.g., '0,0', '0,1', etc.) and corresponding DataFrames as values.
    - subset_dict_data (dict): A dictionary containing descriptive titles for each subset, used to label each subplot.
    - title (str): Title for the entire figure (default is empty).
    - cmap (str): Color map for the heatmaps (default is 'viridis').
    - annotate (bool): Whether to annotate the cells with their values (default is True).
    '''
    
    # Create a 2x2 grid of subplots with a specified figure size
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
   
    # Iterate through the subset dictionary to plot each subset
    for axses, data in subset_dict.items():
        # Extract the relevant columns (TUPOR, SESY, ASER, ACR) for the heatmap
        df = data[['TUPOR', 'SESY', 'ASER']]
        
        # Set the index of the dataframe to the 'name' attribute of the data
        df.index = data.name.tolist()  # Using names as index
        
        # Determine the position of the subplot based on the key (e.g., '0,0', '0,1')
        i = int(axses.split(',')[0])  # Row index
        j = int(axses.split(',')[1])  # Column index
        ax = axes[i, j]

        # Plot the heatmap for the current subset on the specified subplot
        sns.heatmap(df, annot=annotate, cmap=cmap, ax=ax)

        # Modify the y-axis labels for better readability by inserting line breaks
        new_labels = [label.get_text().replace('_epsilon', '\n epsilon').replace('_mut_r', '\n mut_r').replace('addcarbon', 'AddCarbon') for label in ax.get_yticklabels()]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=11)
        if ph4:
            ax.set_xticklabels(labels=['TUPOR_pharm', 'SESY_pharm', 'ASER_pharm'], fontsize=11)
        # Set the title for the current subplot to indicate the subset
        ax.set_title(f"{subset_dict_data[axses]}")
    
    # Set the overall title for the figure
    fig.suptitle(f'{title}', fontsize=16)
    
    # Adjust layout to ensure titles and labels are well placed
    plt.tight_layout()


    if ph4:
        #plt.savefig(f'img_pharm/heat_mapa/{receptor}/{name_save}.svg', format="svg")
        plt.savefig(f'img_pharm/heat_mapa/{receptor}/{name_save}.png', format="png")
    else:
        plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
        plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.png', format="png")
    # Display the heatmap figure
    plt.show()



def plot_heatmaps_with_diff_from_baseline(baseline_df_all, data_dict, type_split, scaf, receptor='', name_save='', poradi = ''):
    """
    Generates heatmaps comparing subsets of data against a baseline.
    The heatmaps highlight the differences from the baseline for values greater than 0.1 or smaller than -0.1.

    Args:
        baseline_df_all (pd.DataFrame): The baseline DataFrame containing reference data.
        data_dict (dict): A dictionary of DataFrames for different subsets (e.g., '_10k', '_100k', '_500k', base).
        type_split (str): The type of data split (e.g., 'dis', 'sim').
        scaf (str): The scaffold type (e.g., 'csk', 'murcko').

    Returns:
        None
    """

    num_subsets = len(data_dict)
    fig, axes = plt.subplots(1, num_subsets, figsize=(12 * num_subsets, 12))  # 1 row, 4 columns for subsets

    baseline_df = baseline_df_all[['TUPOR', 'SESY', 'ASER']]  
    baseline_df.index = baseline_df_all.name.tolist()

    for k, (subset, df) in enumerate(data_dict.items()):
        normalized_df = df[['TUPOR', 'SESY', 'ASER']]
        normalized_df.index = df.name.tolist()

        baseline_df = baseline_df.set_index(normalized_df.index)  # Ensure same order

        ax = axes[k]

        if subset == '':
            diff_df = baseline_df  
        else:
            diff_df = normalized_df - baseline_df  # Compute difference without absolute value
            mask = (diff_df.abs() <= 0.1)  # Mask values in the range [-0.1, 0.1]
            diff_df[mask] = np.nan  # Hide values close to zero
        
        sns.heatmap(diff_df, annot=True, cmap='coolwarm', cbar_kws={'label': 'Difference from Baseline'}, ax=ax, annot_kws={"size": 30})
        ax.figure.axes[-1].yaxis.label.set_size(20)

        if subset == '':
            subset = 'Full OS'
        elif subset == '_62.5k':
            subset = '62,500'
        else:
            subset = subset.replace('_', '').replace('k', ',000')

        if subset == 'Full OS':
            ax.set_title(f"{subset}",  fontsize=35, wrap=True)
        else:
            ax.set_title(f"{subset} subset",  fontsize=35, wrap=True)

        new_labels = [label.get_text().replace('_epsilon', '\n epsilon').replace('_mut_r', '\n mut_r').replace('addcarbon', 'AddCarbon') for label in ax.get_yticklabels()]
        new_labels = [label.replace('_62.5k', '').replace('_125k', '').replace('_250k', '').replace('_500k', '') for label in new_labels]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=30)
        ax.set_xticklabels(ax.get_xticklabels(), ha="center", fontsize=30)
        ax.set_facecolor('white')

    if scaf == 'csk':
        scaf_str = 'CSK'
    else:
        scaf_str = scaf

    fig.text(
    0.005, 0.97, poradi,
    ha='left', va='top',
    fontsize=40
    )

    #fig.suptitle(f'Heatmaps with Differences from Baseline by more than ±0.1 for {scaf_str} scaffolds and {type_split} split for {receptor.replace("_", " ")}', fontsize=40)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.png', format="png")
    plt.show()