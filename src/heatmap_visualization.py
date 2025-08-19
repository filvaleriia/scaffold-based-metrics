import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.colors import ListedColormap



def preprocesing(type_cluster, type_scaffold, generators_name_list, receptor):
    '''
    Function for connection all data set with normalization
    '''
    # Define path to data

    link = f"data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}"

    link_mean = [f"{link}/{generator}/{generator}_mean_{type_scaffold}_{type_cluster}.csv" for generator in generators_name_list]
    
    # Load data
    df_list = [pd.read_csv(f) for f in link_mean]
    df = pd.concat(df_list, axis=0, ignore_index=True)

    scaler = MinMaxScaler()
    numeric_columns = df.select_dtypes(include=['number']).columns  # Select only numeric columns
    df[numeric_columns] = scaler.fit_transform(df[numeric_columns])  # Apply normalization
    df["name"] = df["name"].str.replace("_mean", "", regex=False)

    
    return df



def preprocesing_org(type_cluster, type_scaffold, generators_name_list, receptor):
    '''
    Function for connection all data set without normalization
    '''
    # Define path to data

    link = f"data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}"

    link_mean = [f"{link}/{generator}/{generator}_mean_{type_scaffold}_{type_cluster}.csv" for generator in generators_name_list]
    
    # Load data
    df_list = [pd.read_csv(f) for f in link_mean]
    df = pd.concat(df_list, axis=0, ignore_index=True)

    df["name"] = df["name"].str.replace("_mean", "", regex=False)

    
    return df



def plot_heatmap(type_cluster, type_scaffold, generators_name_list, receptor, title='', name_save='', cmap='viridis', annotate=True, using_norm_values = True):
    ''' 
    Plots a single heatmap for the given data split.
    
    Args:
    - data (pd.DataFrame): The input data containing the columns to be visualized.
    - title (str): The title of the heatmap (default is empty).
    - cmap (str): The color map to be used for the heatmap (default is 'viridis').
    - annotate (bool): Whether to annotate the cells with their values (default is True).
    - using_norm_values: whether to use normalized values
    '''

    # Extract relevant columns (TUPOR, SESY, ASER, ASR) for visualization
    if using_norm_values:
        data = preprocesing(type_cluster, type_scaffold, generators_name_list, receptor)
    else:
        data = preprocesing_org(type_cluster, type_scaffold, generators_name_list, receptor)

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

    plt.savefig(f'img/heat_map/{receptor}/{name_save}.svg', format="svg")
    plt.savefig(f'img/heat_map/{receptor}/{name_save}.png', format="png")
    # Display the heatmap
    plt.show()



def plot_all_subsets(subset_dict, title='', receptor = '', name_save = '', cmap='viridis', annotate=True, numering = ''):
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
    0.005, 0.97, numering,
    ha='left', va='top',
    fontsize=40
    )
    
    # Set the overall title for the figure
    fig.suptitle(f'{title}', fontsize=40)
    
    # Adjust layout to ensure titles and labels are well placed
    plt.tight_layout()
    plt.savefig(f'img/heat_map/{receptor}/{name_save}.svg', format="svg")
    plt.savefig(f'img/heat_map/{receptor}/{name_save}.png', format="png")
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

        # Set the title for the current subplot to indicate the subset
        ax.set_title(f"{subset_dict_data[axses]}")
    
    # Set the overall title for the figure
    fig.suptitle(f'{title}', fontsize=16)
    
    # Adjust layout to ensure titles and labels are well placed
    plt.tight_layout()



    plt.savefig(f'img/heat_map/{receptor}/{name_save}.svg', format="svg")
    plt.savefig(f'img/heat_map/{receptor}/{name_save}.png', format="png")
    # Display the heatmap figure
    plt.show()



def plot_heatmaps_with_diff_from_baseline(baseline_df_all, data_dict, type_split, scaf, receptor='', name_save='', numering = ''):
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
    0.005, 0.97, numering,
    ha='left', va='top',
    fontsize=40
    )

    #fig.suptitle(f'Heatmaps with Differences from Baseline by more than ±0.1 for {scaf_str} scaffolds and {type_split} split for {receptor.replace("_", " ")}', fontsize=40)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f'img/heat_map/{receptor}/{name_save}.svg', format="svg")
    plt.savefig(f'img/heat_map/{receptor}/{name_save}.png', format="png")
    plt.show()



def plot_combined_heatmap(generators, receptors, scaffolds, splits, metrics, cmap="viridis", title=None, save_name="heamps",  using_norm_values=False):
    """
    Create and save combined heatmap for given generators, receptors, scaffolds, and metrics.

    Parameters
    ----------
    hv : object
        Your preprocessing module with hv.preprocesing function
    generators : list
        List of generator names
    receptors : list
        List of receptors
    scaffolds : list
        List of scaffold types
    splits : list
        List of split types (e.g. ['dis', 'sim'])
    metrics : list
        List of metric names
    cmap : str
        Colormap for heatmap (default "viridis")
    title : str
        Title for the heatmap
    save_name : str
        Base name for saving the figure (no extension)
    - using_norm_values:  bool
        whether to use normalized valuess
    """

    # Build dataframe with all values
    data = []
    for gen in generators:
        for receptor in receptors:
            for type_scaffold in scaffolds:
                for type_cluster in splits:
                    if using_norm_values:
                        df = preprocesing(type_cluster, type_scaffold, generators, receptor)
                    else:
                        df = preprocesing_org(type_cluster, type_scaffold, generators, receptor)
                    for met in metrics:
                        value = df[df.name.str.startswith(gen)][met].iloc[0]
                        data.append([gen, receptor, type_scaffold, type_cluster, met, value])

    df = pd.DataFrame(data, columns=['Generator', 'Receptor', 'Scaffold', 'Split', 'Metric', 'Value'])
    
    # Build heatmap data
    heatmap_data = []
    for generator in df.Generator.unique():
        df_sub = df[df.Generator == generator]
        matrix = []
        for receptor in df_sub.Receptor.unique():
            df_sub_r = df_sub[df_sub.Receptor == receptor]
            for scaffold in df_sub_r.Scaffold.unique():
                df_sub_r_s = df_sub_r[df_sub_r.Scaffold == scaffold]
                for split in splits:
                    df_sub_r_s_s = df_sub_r_s[df_sub_r_s.Split == split]
                    for met in metrics:
                        matrix.append(df_sub_r_s_s[df_sub_r_s_s.Metric == met]['Value'].values[0])
        heatmap_data.append(matrix)

    heatmap_data = np.array(heatmap_data)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(20, 10))
    sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap=cmap, ax=ax,
                cbar_kws={'label': 'Metric Value'}, annot_kws={"size": 12})

    if title:
        ax.set_title(title, fontsize=20)
    ax.set_ylabel('Generators', fontsize=15)

    new_labels = [label.replace('_epsilon', '\n epsilon')
                        .replace('_mut_r', '\n mut_r')
                        .replace('addcarbon', 'AddCarbon')
                  for label in generators]
    ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=15)
    ax.tick_params(axis='y', pad=15)

    # Add rectangles to highlight blocks
    for i in range(len(heatmap_data)):
        for k in [0,3,6,9,12,15,18,21]:
            ax.add_patch(plt.Rectangle((k, i), 3, 1, fill=False, edgecolor='black', lw=1))
    ax.vlines(x=12, ymin=0, ymax=len(generators) * 2, colors='red', linewidth=5)  

    custom_xticklabels = ["CSK-DIS",  "CSK-SIM", "MURCKO-DIS",  "MURCKO-SIM",  "CSK-DIS",  "CSK-SIM", "MURCKO-DIS",  "MURCKO-SIM"]

    x = -1.5
    tick_positions = []
    for i in range(len(custom_xticklabels)):
        x += 3
        tick_positions.append(x)


    ax.set_xticks(tick_positions)
    ax.set_xticklabels(custom_xticklabels, rotation=45, ha="right", fontsize=15)

    ax.text(2, 11, "Glucocorticoid receptor", fontsize=17, color='black')
    ax.text(17, 11, "Leukocyte elastase", fontsize=17, color='black')

    plt.tight_layout()
    plt.savefig(f'img/heat_map/{save_name}.svg', format="svg")
    plt.savefig(f'img/heat_map/{save_name}.png', format="png")
    plt.show()



def make_cmap_to_white(base_hex_color):
    # Convert a base hex color to RGB
    base_rgb = to_rgb(base_hex_color)
    # Define white color in RGB
    white_rgb = to_rgb('#f0f0f0')
    # Create a gradient from white to the base color
    colors = [white_rgb, base_rgb]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    return cmap



def plot_combined_heatmap_variable_cmaps(
        generators, receptors, scaffolds, splits,
        metrics=['TUPOR', 'SESY', 'ASER'], 
        title=None, save_name=None, using_norm_values=False):
    """
    Plot combined heatmaps for multiple generators, receptors, scaffolds, and metrics.
    
    Parameters:
    - generators: list of generator names
    - receptors: list of receptor names
    - scaffolds: list of scaffold types
    - splits: list of data splits
    - metrics: list of metrics to plot
    - title: figure title
    - save_name: name for saving the figure
    - using_norm_values: whether to use normalized values
    """
    
    # Define base colors for each metric
    metric_base_colors = {
        'TUPOR': "#e97b32",
        'SESY': "#97C2F0",  
        'ASER': "#71ad48"
    }

    # Collect all metric values into a single DataFrame
    data = []
    for gen in generators:
        for receptor in receptors:
            for type_scaffold in scaffolds:
                for type_cluster in splits:
                    if using_norm_values:
                        df = preprocesing(type_cluster, type_scaffold, generators, receptor)
                    else:
                        df = preprocesing_org(type_cluster, type_scaffold, generators, receptor)
                    for met in metrics:
                        value = df[df.name.str.startswith(gen)][met].iloc[0]
                        data.append([gen, receptor, type_scaffold, type_cluster, met, value])

    df = pd.DataFrame(data, columns=['Generator', 'Receptor', 'Scaffold', 'Split', 'Metric', 'Value'])

    # Create subplots dynamically based on number of receptors
    nrows = len(receptors)
    ncols = len(metrics)
    fig_width = max(5*ncols, 8)
    fig_height = max(5*nrows, 3*nrows + 1)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width, fig_height))

    # Ensure axes is always 2D array for consistent indexing
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])

    # Loop through metrics and receptors
    for col_idx, metric in enumerate(metrics):
        metric_df = df[df['Metric'] == metric].copy()

        for row_idx, receptor in enumerate(receptors):
            ax = axes[row_idx, col_idx]

            # Prepare heatmap data
            heatmap_data = []
            for gen in generators:
                row = []
                for scaffold in scaffolds:
                    for split in splits:
                        mask = (
                            (metric_df['Generator'] == gen) &
                            (metric_df['Receptor'] == receptor) &
                            (metric_df['Scaffold'] == scaffold) &
                            (metric_df['Split'] == split)
                        )
                        value = metric_df[mask]['Value'].values[0]
                        row.append(value)
                heatmap_data.append(row)

            heatmap_array = np.array(heatmap_data)
            cmap_custom = make_cmap_to_white(metric_base_colors[metric])

            # Plot heatmap with custom colors and annotation
            if using_norm_values:
                sns.heatmap(
                    heatmap_array, annot=True, cmap=cmap_custom, ax=ax,
                    cbar_kws={'label': metric}, annot_kws={"size": 11, "color": "black"},
                    vmin=metric_df['Value'].min(), vmax=metric_df['Value'].max()
                )
            else:
                sns.heatmap(
                    heatmap_array, annot=True, fmt=".4f", cmap=cmap_custom, ax=ax,
                    cbar_kws={'label': metric}, annot_kws={"size": 11, "color": "black"},
                    vmin=metric_df['Value'].min(), vmax=metric_df['Value'].max()
                )

            # Set titles and y-axis labels
            if row_idx == 0:
                ax.set_title(metric, fontsize=16)
            if col_idx == 0:
                ax.set_ylabel(receptor.replace("_", " "), fontsize=14)
                ax.set_yticks(np.arange(len(generators)) + 0.5)
                new_labels = [label.replace('_epsilon', '\n epsilon')
                                   .replace('_mut_r', '\n mut_r')
                                   .replace('addcarbon', 'AddCarbon')
                              for label in generators]
                ax.set_yticklabels(new_labels, rotation=0, fontsize=10)
            else:
                ax.set_ylabel("")
                ax.set_yticks([])
                ax.set_yticklabels([])

            # X-axis labels for bottom row only
            if row_idx == nrows-1:
                xticklabels = [f"{sc}-{split}" for sc in scaffolds for split in splits]
                ax.set_xticks(np.arange(len(xticklabels)) + 0.5)
                ax.set_xticklabels(xticklabels, rotation=45, ha="right", fontsize=10)
            else:
                ax.set_xticks(np.arange(len(scaffolds)*len(splits)) + 0.5)
                ax.set_xticklabels([])

    # Set overall figure title
    if title:
        fig.suptitle(title, fontsize=18)

    # Adjust layout and save
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    if save_name:
        plt.savefig(f'img/heat_map/{save_name}.svg', format="svg")
        plt.savefig(f'img/heat_map/{save_name}.png', format="png")

    plt.show()



def plot_combined_heatmap_with_single_column_for_each_metric(
        generators, receptors, scaffolds, splits,
        metrics=['TUPOR', 'SESY', 'ASER'], 
        title=None, save_name=None, using_norm_values=False,
        inter_metric_wspace=0.15,   # Larger spacing between metrics
        intra_metric_wspace=0.05    # Smaller spacing within a metric block
    ):

    """
    Plot combined heatmaps for multiple generators, receptors, scaffolds, and metrics.
    
    Parameters:
    - generators: list of generator names
    - receptors: list of receptor names
    - scaffolds: list of scaffold types
    - splits: list of data splits
    - metrics: list of metrics to plot
    - title: figure title
    - save_name: name for saving the figure
    - using_norm_values: whether to use normalized values
    - inter_metric_wspace: space between different metrics
    - intra_metric_wspace: space between individuals scaffolds and split in one Metric
    """

    # Base colors for each metric
    metric_base_colors = {
        'TUPOR': "#e97b32",
        'SESY': "#97C2F0",
        'ASER': "#71ad48"
    }

    # --- Build a DataFrame with all values ---
    data = []
    for gen in generators:
        for receptor in receptors:
            for type_scaffold in scaffolds:
                for type_cluster in splits:
                    # Select preprocessed data, normalized or original
                    df = (preprocesing_org if not using_norm_values else preprocesing)(
                        type_cluster, type_scaffold, generators, receptor
                    )
                    # Extract values for each metric
                    for met in metrics:
                        value = df[df.name.str.startswith(gen)][met].iloc[0]
                        data.append([gen, receptor, type_scaffold, type_cluster, met, value])

    # Convert the collected data into a pandas DataFrame
    df = pd.DataFrame(
        data, columns=['Generator', 'Receptor', 'Scaffold', 'Split', 'Metric', 'Value']
    )

    nrows = len(receptors)
    nmetrics = len(metrics)

    # Figure size: wider figure depending on number of metrics + spacing
    fig_width = 1.7 * 4 * nmetrics + 2
    fig_height = 6 * nrows
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Outer grid for arranging metrics per receptor row
    outer_gs = fig.add_gridspec(
        nrows=nrows, ncols=nmetrics,
        wspace=inter_metric_wspace, hspace=0.1
    )

    for met_idx, metric in enumerate(metrics):
        # Filter data for the current metric
        metric_df = df[df['Metric'] == metric].copy()
        cmap_custom = make_cmap_to_white(metric_base_colors[metric])

        for rec_idx, receptor in enumerate(receptors):
            # Inner sub-grid for each metric and receptor
            inner = outer_gs[rec_idx, met_idx].subgridspec(
                nrows=1, ncols=4, wspace=intra_metric_wspace, hspace=0.0
            )

            group_axes = []  # store 4 axes for top labels
            # Plot 4 heatmaps: csk-[dis,sim], murcko-[dis,sim]
            for sc_idx, scaffold_type in enumerate(["csk", "murcko"]):
                # Extract block data for this scaffold type
                block_df = metric_df[
                    (metric_df['Receptor'] == receptor) &
                    (metric_df['Scaffold'] == scaffold_type)
                ]
                
                for split_idx, split in enumerate(["dis", "sim"]):
                    col = sc_idx * 2 + split_idx
                    ax = fig.add_subplot(inner[0, col])
                    group_axes.append(ax)

                    # Prepare heatmap array
                    sub_df = block_df[block_df['Split'] == split].copy()
                    sub_df = sub_df.set_index('Generator').reindex(generators)
                    heatmap_array = sub_df['Value'].to_numpy().reshape(-1, 1)

                    vmin = heatmap_array.min()
                    vmax = heatmap_array.max()

                    # Show colorbar only on the last subplot of the group
                    show_colorbar = (sc_idx == 1 and split_idx == 1)

                    if using_norm_values:
                        sns.heatmap(
                            heatmap_array,
                            annot=True,
                            cmap=cmap_custom, ax=ax,
                            cbar=show_colorbar,
                            cbar_kws={'label': metric} if show_colorbar else None,
                            annot_kws={"size": 13, "color": "black"},
                            vmin=vmin, vmax=vmax
                        )
                    else:
                        sns.heatmap(
                            heatmap_array,
                            annot=True, fmt=".4f",
                            cmap=cmap_custom, ax=ax,
                            cbar=show_colorbar,
                            cbar_kws={'label': metric} if show_colorbar else None,
                            annot_kws={"size": 13, "color": "black"},
                            vmin=vmin, vmax=vmax
                        )
                    ax.set_aspect("auto")

                    # X-axis = split
                    ax.set_xticks([0.5])
                    ax.set_xticklabels([split], rotation=0, ha="center", fontsize=12)

                    # Y-axis only for the first metric and first scaffold (leftmost)
                    if met_idx == 0 and sc_idx == 0 and split_idx == 0:
                        ax.set_ylabel(receptor.replace('_', ' '), fontsize=13)
                        ax.set_yticks(np.arange(len(generators)) + 0.5)
                        new_labels = [g.replace('_epsilon', '\n epsilon')
                                        .replace('_mut_r', '\n mut_r')
                                        .replace('addcarbon', 'AddCarbon')
                                      for g in generators]
                        ax.set_yticklabels(new_labels, rotation=0, fontsize=13)
                    else:
                        ax.set_ylabel("")
                        ax.set_yticks([])
                        ax.set_yticklabels([])

            # Top labels above each scaffold pair (only for top receptor row)
            if rec_idx == 0:
                # CSK (columns 0 and 1)
                p0 = group_axes[0].get_position()
                p1 = group_axes[1].get_position()
                x_mid_csk = (p0.x0 + p1.x1) / 2
                y_top_csk = max(p0.y1, p1.y1) + 0.012
                fig.text(x_mid_csk, y_top_csk, f"{metric} - CSK", ha="center", va="bottom", fontsize=13)

                # MURCKO (columns 2 and 3)
                p2 = group_axes[2].get_position()
                p3 = group_axes[3].get_position()
                x_mid_mur = (p2.x0 + p3.x1) / 2
                y_top_mur = max(p2.y1, p3.y1) + 0.012
                fig.text(x_mid_mur, y_top_mur, f"{metric} - MURCKO", ha="center", va="bottom", fontsize=13)

    # Add a global title if specified
    if title:
        fig.suptitle(title, fontsize=14, y=0.995)

    # Adjust layout and save the figure
    plt.tight_layout(rect=[0, 0, 1, 0.965])
    if save_name:
        plt.savefig(f'img/heat_map/{save_name}.svg', format="svg", bbox_inches='tight')
        plt.savefig(f'img/heat_map/{save_name}.png', format="png", dpi=200, bbox_inches='tight')
    plt.show()


