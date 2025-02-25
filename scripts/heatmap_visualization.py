import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd



def plot_heatmap(data, title='', name_save='',receptor = '', cmap='viridis', annotate=True):
    ''' 
    Plots a single heatmap for the given data split.
    
    Args:
    - data (pd.DataFrame): The input data containing the columns to be visualized.
    - title (str): The title of the heatmap (default is empty).
    - cmap (str): The color map to be used for the heatmap (default is 'viridis').
    - annotate (bool): Whether to annotate the cells with their values (default is True).
    '''

    # Extract relevant columns (TUPOR, SESY, ASER, ASR) for visualization
    df = data[['TUPOR', 'SESY', 'ASER', 'ACR']]
    
    # Set the index of the dataframe to the 'name' attribute of the data
    df.index = data.name.tolist()

    # Create a figure for the heatmap with a specific size
    plt.figure(figsize=(10, 6))
    
    # Plot the heatmap using seaborn with optional annotations and custom color map
    sns.heatmap(df, annot=annotate, cmap=cmap)
    
    # Set the title for the heatmap
    plt.title(title)
    plt.tight_layout()
    # Save the plot as an SVG file
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
    # Display the heatmap
    plt.show()



def plot_all_subsets(subset_dict, title='', receptor = '', name_save = '', cmap='viridis', annotate=True):
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
    fig, axes = plt.subplots(1, num_subsets, figsize=(num_subsets * 6, 6))
    
    # If there is only one subset, axes will not be iterable, so convert it to a list
    if num_subsets == 1:
        axes = [axes]
    
    # Iterate over each axis (subplot) and the corresponding subset data
    for ax, (subset_name, data) in zip(axes, subset_dict.items()):
        # Extract the relevant columns for the heatmap (TUPOR, SESY, ASER, ACR)
        df = data[['TUPOR', 'SESY', 'ASER', 'ACR']]
        
        # Set the index of the dataframe to the 'name' attribute of the data
        df.index = data.name.tolist()  # Using names as index
        
        # Plot the heatmap for the current subset
        sns.heatmap(df, annot=annotate, cmap=cmap, ax=ax)
        
        # If the subset name is empty, label it as 'base'
        if subset_name == '':
            subset_name = 'base'
        
        # Modify the y-axis labels for better readability by replacing certain substrings
        new_labels = [label.get_text().replace('_epsilon', '\n epsilon').replace('_mean', '\n mean') for label in ax.get_yticklabels()]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=11)
        
        # Set the title for the current subplot to indicate the subset name
        ax.set_title(f"{subset_name} subset")
    
    # Set the overall title for the figure
    fig.suptitle(f'{title}', fontsize=16)
    
    # Adjust layout to ensure titles and labels are well placed
    plt.tight_layout()
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
    # Display the heatmap figure
    plt.show()



def plot_heatmap_base(subset_dict, subset_dict_data, title='', receptor = '', name_save = '', cmap='viridis', annotate=True):
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
        df = data[['TUPOR', 'SESY', 'ASER', 'ACR']]
        
        # Set the index of the dataframe to the 'name' attribute of the data
        df.index = data.name.tolist()  # Using names as index
        
        # Determine the position of the subplot based on the key (e.g., '0,0', '0,1')
        i = int(axses.split(',')[0])  # Row index
        j = int(axses.split(',')[1])  # Column index
        ax = axes[i, j]

        # Plot the heatmap for the current subset on the specified subplot
        sns.heatmap(df, annot=annotate, cmap=cmap, ax=ax)

        # Modify the y-axis labels for better readability by inserting line breaks
        new_labels = [label.get_text().replace('_epsilon', '\n epsilon').replace('_mean', '\n mean') for label in ax.get_yticklabels()]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=11)
        
        # Set the title for the current subplot to indicate the subset
        ax.set_title(f"{subset_dict_data[axses]}")
    
    # Set the overall title for the figure
    fig.suptitle(f'{title}', fontsize=16)
    
    # Adjust layout to ensure titles and labels are well placed
    plt.tight_layout()
    # Save the plot as an SVG file
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
    # Display the heatmap figure
    plt.show()



def plot_heatmaps_from_baseline(base_line, data_dict, type_split, scaf, receptor = '', name_save = ''):
    """
    Plots heatmaps comparing different subsets to a baseline.
    
    Args:
    - base_line (pd.DataFrame): DataFrame containing the baseline data to compare against.
    - data_dict (dict): Dictionary with keys representing subset names (e.g., '_10k', '_100k', '_500k', '') 
                        and values representing corresponding DataFrames.
    """

    # Create figure for subplots (1 row, 4 columns for subsets)
    fig, axes = plt.subplots(1, 4, figsize=(30, 8))  # 1 row, 4 columns for subsets
    # Set baseline for comparison
    baseline_df = base_line[['TUPOR', 'SESY', 'ASER', 'ACR']].copy()
    baseline_df.index = base_line.name.tolist()

    # Iterate over subsets
    for k, (subset, df) in enumerate(data_dict.items()):
        # Load data for current subset
        normalized_df = df[['TUPOR', 'SESY', 'ASER', 'ACR']]
        normalized_df.index = df.name.tolist()
        # Ensure that baseline and current subset have the same index order and columns
        baseline_df = baseline_df.set_index(normalized_df.index)
        ax = axes[k]  # Select subplot
        # For baseline, display all values; for others, apply mask
        if subset == '':
            mask = np.zeros_like(normalized_df, dtype=bool)  # No masking (display all)
        else:
            mask = np.abs(normalized_df - baseline_df) < 0.1  # Mask values close to baseline
        # Plot heatmap with mask
        sns.heatmap(normalized_df, annot=True, cmap='viridis', cbar_kws={'label': 'Normalized Value'}, ax=ax, mask=mask)
        # Set title for the subset
        ax.set_title(f'{scaf} - {subset if subset else "BASELINE"}', fontsize=10, wrap=True)
        # Wrap generator names for readability
        new_labels = [label.get_text().replace('_epsilon', '\n epsilon').replace('_mean', '\n mean') for label in ax.get_yticklabels()]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=11)
    # Add overall title
    fig.suptitle(f'Heatmaps with Differences from Baseline > 0.1 ({scaf}, {type_split}) Original values for {receptor}', fontsize=16)
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])  
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
    plt.show()



def plot_heatmaps_with_diff_from_baseline(baseline_df_all, data_dict, type_split, scaf, receptor='', name_save=''):
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
    fig, axes = plt.subplots(1, 4, figsize=(4 * 6, 6))  # 1 row, 4 columns for subsets

    baseline_df = baseline_df_all[['TUPOR', 'SESY', 'ASER', 'ACR']]  
    baseline_df.index = baseline_df_all.name.tolist()

    for k, (subset, df) in enumerate(data_dict.items()):
        normalized_df = df[['TUPOR', 'SESY', 'ASER', 'ACR']]
        normalized_df.index = df.name.tolist()

        baseline_df = baseline_df.set_index(normalized_df.index)  # Ensure same order

        ax = axes[k]

        if subset == '':
            diff_df = baseline_df  
        else:
            diff_df = normalized_df - baseline_df  # Compute difference without absolute value
            mask = (diff_df.abs() <= 0.1)  # Mask values in the range [-0.1, 0.1]
            diff_df[mask] = np.nan  # Hide values close to zero

        sns.heatmap(diff_df, annot=True, cmap='coolwarm', cbar_kws={'label': 'Difference from Baseline'}, ax=ax)

        ax.set_title(f'{scaf} - {subset if subset else "BASELINE"}', fontsize=11, wrap=True)

        new_labels = [label.get_text().replace('_epsilon', '\n epsilon').replace('_mean', '\n mean') for label in ax.get_yticklabels()]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=11)

    fig.suptitle(f'Heatmaps with Differences from Baseline > 0.1 or < -0.1 ({scaf}, {type_split})', fontsize=16)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
    plt.show()