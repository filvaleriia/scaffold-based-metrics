import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
from matplotlib.gridspec import GridSpec


def preprocesing(type_cluster, type_scaffold, generators_name_list, receptor, data_folder):
    '''
    Function for connection all data set with normalization
    '''
    # Define path to data

    link = f"{data_folder}data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}"

    link_mean = [f"{link}/{generator}/{generator}_mean_{type_scaffold}_{type_cluster}.csv" for generator in generators_name_list]
    
    # Load data
    df_list = [pd.read_csv(f) for f in link_mean]
    df = pd.concat(df_list, axis=0, ignore_index=True)

    scaler = MinMaxScaler()
    numeric_columns = df.select_dtypes(include=['number']).columns  # Select only numeric columns
    df[numeric_columns] = scaler.fit_transform(df[numeric_columns])  # Apply normalization
    df["name"] = df["name"].str.replace("_mean", "", regex=False)

    
    return df



def preprocesing_org(type_cluster, type_scaffold, generators_name_list, receptor, data_folder):
    '''
    Function for connection all data set without normalization
    '''
    # Define path to data

    link = f"{data_folder}data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}"
    link_mean = [f"{link}/{generator}/{generator}_mean_{type_scaffold}_{type_cluster}.csv" for generator in generators_name_list]
    
    # Load data
    df_list = [pd.read_csv(f) for f in link_mean]
    df = pd.concat(df_list, axis=0, ignore_index=True)

    df["name"] = df["name"].str.replace("_mean", "", regex=False)

    
    return df



def plot_heatmap(type_cluster, type_scaffold, generators_name_list, receptor, title='', name_save='', cmap='viridis', annotate=True, using_norm_values = True, data_folder = '', save_folder = ''):
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
        data = preprocesing(type_cluster, type_scaffold, generators_name_list, receptor, data_folder)
    else:
        data = preprocesing_org(type_cluster, type_scaffold, generators_name_list, receptor, data_folder)

    df = data[['TUPOR', 'SESY', 'ASER']]
    # Set the index of the dataframe to the 'name' attribute of the data
    df.index = data.name.tolist()

    # Create a figure for the heatmap with a specific size
    plt.figure(figsize=(10, 0.7 * len(generators_name_list)))
    
    # Plot the heatmap using seaborn with optional annotations and custom color map
    sns.heatmap(df, annot=annotate, cmap=cmap,   annot_kws={"size": 17})
    
    # Set the title for the heatmap
    plt.title(title, fontsize=17 , pad=10)

    new_labels = [label.replace('_epsilon', '\n epsilon').replace('_mut_r', '\n mut_r').replace('addcarbon', 'AddCarbon') for label in df.index]
    #new_labels = [label.replace('_62.5k', '').replace('_125k', '').replace('_250k', '').replace('_500k', '') for label in new_labels]
    plt.xticks(fontsize=17)
    plt.yticks(ticks=np.arange(len(df.index)) + 0.5,labels=new_labels,fontsize=17)
    plt.tight_layout()
    # Save the plot as an SVG file
    if save_folder:
        plt.savefig(f'{save_folder}/{name_save}.svg', format="svg")
        plt.savefig(f'{save_folder}/{name_save}.png', format="png")
    else:
        plt.savefig(f'img/heat_map/{receptor}/{name_save}.svg', format="svg")
        plt.savefig(f'img/heat_map/{receptor}/{name_save}.png', format="png")
    # Display the heatmap
    plt.show()



def plot_all_subsets(subset_dict, title='', receptor='', name_save='', cmap='viridis', annotate=True, numering='', save_folder=''):
    '''
    Plots heatmaps for multiple subsets in a single figure in paper style,
    with proper spacing and colorbar handling as in plot_heatmaps_with_diff_from_baseline.
    '''

    num_subsets = len(subset_dict)
    num_gen = len(subset_dict[next(iter(subset_dict))])

    # Dynamically determine columns and rows
    max_cols = 5
    cols = min(num_subsets, max_cols)
    rows = math.ceil(num_subsets / cols)

    # Adjust layout to avoid sparse last row
    if num_subsets > 1:
        cols_last_row = num_subsets % cols
        if cols_last_row != 0 and cols_last_row < math.ceil(cols / 2):
            cols = math.ceil(num_subsets / 2)
            rows = math.ceil(num_subsets / cols)

    # Figure size in paper style
    fig_width = cols * 12
    fig_height = rows * num_gen * 1.3
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(rows, cols, figure=fig)

    axes = []
    for i in range(num_subsets):
        row = i // cols
        col = i % cols
        ax = fig.add_subplot(gs[row, col])
        axes.append(ax)

    # Plot each heatmap
    for ax, (subset_name, data) in zip(axes, subset_dict.items()):
        df = data[['TUPOR', 'SESY', 'ASER']]
        df.index = data.name.tolist()

        hm = sns.heatmap(
            df,
            annot=annotate,
            cmap=cmap,
            ax=ax,
            annot_kws={"size": 30},  # paper style
            cbar=False
        )

        divider = make_axes_locatable(ax)
        cax = divider.append_axes(
            "right",
            size="2.0%",
            pad=0.15
        )
        cbar = plt.colorbar(hm.collections[0], cax=cax)
        cbar.ax.tick_params(labelsize=15)
        cbar.outline.set_visible(False)
        ax.figure.axes[-1].yaxis.label.set_size(15)

        # Format subset name
        if subset_name == '':
            subset_name_disp = 'Full OS'
        elif subset_name == '_62.5k':
            subset_name_disp = '62,500'
        else:
            subset_name_disp = subset_name.replace('_', '').replace('k', ',000')

        ax.set_title(f"{subset_name_disp}" if subset_name_disp == 'Full OS' else f"{subset_name_disp} subset",
                     fontsize=35, wrap=True, pad=12)

        # Format y-axis labels
        new_labels = [
            label.get_text()
            .replace('_epsilon', '\n epsilon')
            .replace('_mut_r', '\n mut_r')
            .replace('addcarbon', 'AddCarbon')
            for label in ax.get_yticklabels()
        ]
        new_labels = [
            label.replace('_62.5k', '').replace('_125k', '')
                 .replace('_250k', '').replace('_500k', '')
                 .replace('_10k', '')
            for label in new_labels
        ]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=30)
        ax.set_xticklabels(ax.get_xticklabels(), ha="center", fontsize=30)
        ax.tick_params(axis='y', pad=12)
        ax.set_facecolor('white')

    # Hide unused subplot slots
    total_slots = rows * cols
    for j in range(len(axes), total_slots):
        row = j // cols
        col = j % cols
        ax_empty = fig.add_subplot(gs[row, col])
        ax_empty.axis('off')

    # Posun druhého řádku pro větší vertikální mezeru
    if rows == 2:
        for i, ax in enumerate(axes):
            if i >= cols:  # druhý řádek
                pos = ax.get_position()
                ax.set_position([pos.x0, pos.y0 + 0.03, pos.width + 0.01, pos.height])


    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Numbering and global title
    fig.text(0.005, 0.95, numering, ha='left', va='top', fontsize=40)
    if title:
        fig.suptitle(title, fontsize=40)

    # Save
    if save_folder:
        plt.savefig(f'{save_folder}/{name_save}.svg', format="svg")
        plt.savefig(f'{save_folder}/{name_save}.png', format="png", dpi=300, bbox_inches="tight")
        plt.savefig(f'{save_folder}/{name_save}.pdf', bbox_inches='tight')
    else:
        plt.savefig(f'img/heat_map/{receptor}/{name_save}.svg', format="svg")
        plt.savefig(f'img/heat_map/{receptor}/{name_save}.png', format="png", dpi=300, bbox_inches="tight")

    plt.show()





def plot_heatmap_base(subset_dict, subset_dict_data, title='', receptor = '', name_save = '', cmap='viridis', annotate=True, save_folder = ''):
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
    fig, axes = plt.subplots(2, 2, figsize=(14, 1.3 * len((subset_dict[next(iter(subset_dict))]))))
   
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


    if save_folder:
        plt.savefig(f'{save_folder}/{name_save}.svg', format="svg")
        plt.savefig(f'{save_folder}/{name_save}.png', format="png")
    else:
        plt.savefig(f'img/heat_map/{receptor}/{name_save}.svg', format="svg")
        plt.savefig(f'img/heat_map/{receptor}/{name_save}.png', format="png")
    # Display the heatmap figure
    plt.show()



def plot_heatmaps_with_diff_from_baseline(baseline_df_all, data_dict, type_split, scaf, receptor='', name_save='', numering='', save_folder=''):
    '''
    Plots heatmaps showing differences from the baseline for multiple subsets,
    dynamically arranging subplots based on the number of subsets.
    '''
    num_subsets = len(data_dict)
    num_gen = len(data_dict[next(iter(data_dict))])

    # Dynamically determine number of rows and columns (max 5 columns)
    max_cols = 5
    cols = min(num_subsets, max_cols)
    rows = math.ceil(num_subsets / cols)

    # Adjust layout so that the last row is not too sparse
    if num_subsets > 1:
        cols_last_row = num_subsets % cols
        if cols_last_row != 0 and cols_last_row < math.ceil(cols / 2):
            cols = math.ceil(num_subsets / 2)
            rows = math.ceil(num_subsets / cols)

    # Figure size in paper style
    fig_width = cols * 12
    fig_height = rows * num_gen * 1.3
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(rows, cols, figure=fig)

    axes = []
    for i in range(num_subsets):
        row = i // cols
        col = i % cols
        ax = fig.add_subplot(gs[row, col])
        axes.append(ax)

    # Prepare baseline values
    baseline_df = baseline_df_all[['TUPOR', 'SESY', 'ASER']]
    baseline_df.index = baseline_df_all.name.tolist()

    # Plot differences
    for ax, (subset, df) in zip(axes, data_dict.items()):
        normalized_df = df[['TUPOR', 'SESY', 'ASER']]
        normalized_df.index = df.name.tolist()
        
        # Align indexes with baseline
        baseline_df = baseline_df.set_index(normalized_df.index)

        # Compute differences; empty subset uses direct baseline values
        if subset == '':
            diff_df = baseline_df
        else:
            diff_df = normalized_df - baseline_df
            mask = (diff_df.abs() <= 0.1)
            diff_df[mask] = np.nan

        hm = sns.heatmap(
            diff_df,
            annot=True,
            cmap='coolwarm',
            cbar=False,
            ax=ax,
            annot_kws={"size": 30}  # paper style
        )

        divider = make_axes_locatable(ax)
        cax = divider.append_axes(
            "right",
            size="2.0%",
            pad=0.15
        )
        cbar = plt.colorbar(hm.collections[0], cax=cax)
        cbar.ax.tick_params(labelsize=15)
        cbar.outline.set_visible(False)
        ax.figure.axes[-1].yaxis.label.set_size(15)

        # Format subset name for title
        if subset == '':
            subset_name = 'Full OS'
        elif subset == '_62.5k':
            subset_name = '62,500'
        else:
            subset_name = subset.replace('_', '').replace('k', ',000')

        title_text = f"{subset_name}" if subset_name == 'Full OS' else f"{subset_name} subset"
        ax.set_title(title_text, fontsize=35, wrap=True, pad=12)

        # Format Y-axis labels
        new_labels = [
            label.get_text()
            .replace('_epsilon', '\n epsilon')
            .replace('_mut_r', '\n mut_r')
            .replace('addcarbon', 'AddCarbon')
            for label in ax.get_yticklabels()
        ]
        new_labels = [
            label.replace('_62.5k', '').replace('_125k', '')
                .replace('_250k', '').replace('_500k', '')
                .replace('_10k', '')
            for label in new_labels
        ]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=30)
        ax.tick_params(axis='y', pad=12)
        ax.set_xticklabels(ax.get_xticklabels(), ha="center",  fontsize=30)
        ax.set_facecolor('white')

    # Hide unused subplot slots
    total_slots = rows * cols
    for j in range(len(axes), total_slots):
        row = j // cols
        col = j % cols
        ax_empty = fig.add_subplot(gs[row, col])
        ax_empty.axis('off')


    if rows == 2:
        for i, ax in enumerate(axes):
            if i >= cols: 
                pos = ax.get_position()
                ax.set_position([pos.x0, pos.y0 + 0.03, pos.width + 0.01, pos.height])

    # Add numbering text
    fig.text(0.005, 0.95, numering, ha='left', va='top', fontsize=40)

    # Format scaffold label
    scaf_str = 'CSK' if scaf == 'csk' else scaf
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save
    if save_folder:
        plt.savefig(f'{save_folder}/{name_save}.svg', format="svg")
        plt.savefig(f'{save_folder}/{name_save}.png', format="png", dpi=300, bbox_inches="tight")
        plt.savefig(f'{save_folder}/{name_save}.pdf', bbox_inches='tight') 
    else:
        plt.savefig(f'img/heat_map/{receptor}/{name_save}.svg', format="svg")
        plt.savefig(f'img/heat_map/{receptor}/{name_save}.png', format="png", dpi=300, bbox_inches="tight")

    plt.show()




def plot_combined_heatmap(generators, receptors, scaffolds, splits, metrics, cmap="viridis", title=None, save_name="heamps",  using_norm_values=False, data_folder = '', save_folder = ''):
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
                        df = preprocesing(type_cluster, type_scaffold, generators, receptor, data_folder)
                    else:
                        df = preprocesing_org(type_cluster, type_scaffold, generators, receptor, data_folder)
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

    if save_folder:
        plt.savefig(f'{save_folder}/{save_name}.svg', format="svg")
        plt.savefig(f'{save_folder}/{save_name}.png', format="png")
    else:
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
        title=None, save_name=None, using_norm_values=False, data_folder='', save_folder = ''):
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
                        df = preprocesing(type_cluster, type_scaffold, generators, receptor, data_folder)
                    else:
                        df = preprocesing_org(type_cluster, type_scaffold, generators, receptor, data_folder)
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
        if save_folder:
            plt.savefig(f'{save_folder}/{save_name}.svg', format="svg")
            plt.savefig(f'{save_folder}/{save_name}.png', format="png")
        else:
            plt.savefig(f'img/heat_map/{save_name}.svg', format="svg")
            plt.savefig(f'img/heat_map/{save_name}.png', format="png")

    plt.show()



def plot_combined_heatmap_with_single_column_for_each_metric(
        generators, receptors, scaffolds, splits,
        metrics=['TUPOR', 'SESY', 'ASER'], 
        title=None, save_name=None, using_norm_values=False, data_folder='', save_folder = '',
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
                        type_cluster, type_scaffold, generators, receptor, data_folder
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
    fig_height = 6 * nrows * 1.3

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

                    heatmap_flat = heatmap_array.flatten()
                    max_idx = np.argmax(heatmap_flat)

                    # Build annot_array with same shape as heatmap_array
                    annot_array = []
                    for i, val in enumerate(heatmap_flat):
                        text = f"{val:.4f}" if not using_norm_values else f"{val:.3f}"
                        if i == max_idx:
                            text = r"$\bf{" + text + "}$"  # bold max
                        annot_array.append(text)

                    # Reshape to heatmap_array shape
                    annot_array = np.array(annot_array).reshape(heatmap_array.shape)

                    # Show colorbar only on the last subplot of the group
                    show_colorbar = (sc_idx == 1 and split_idx == 1)
                    if show_colorbar:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                    else:
                        cax = None

                    if using_norm_values:
                        sns.heatmap(
                            heatmap_array,
                            annot=annot_array,
                            fmt="",
                            cmap=cmap_custom, ax=ax,
                            cbar=show_colorbar,
                            cbar_ax=cax,
                            cbar_kws={'label': metric} if show_colorbar else None,
                            annot_kws={"size": 15, "color": "black"},
                            vmin=vmin, vmax=vmax
                        )
                    else:
                        sns.heatmap(
                            heatmap_array,
                            annot=annot_array, 
                            fmt="",
                            cmap=cmap_custom, ax=ax,
                            cbar=show_colorbar,
                            cbar_ax=cax,
                            cbar_kws={'label': metric} if show_colorbar else None,
                            annot_kws={"size": 14, "color": "black"},
                            vmin=vmin, vmax=vmax
                        )
                    ax.set_aspect("auto")

                    # X-axis = split
                    ax.set_xticks([0.5])
                    ax.set_xticklabels([split], rotation=0, ha="center", fontsize=12)

                    # Y-axis only for the first metric and first scaffold (leftmost)
                    if met_idx == 0 and sc_idx == 0 and split_idx == 0:
                        ax.set_ylabel(receptor.replace('_', ' '), fontsize=14)
                        ax.set_yticks(np.arange(len(generators)) + 0.5)
                        new_labels = [g.replace('_epsilon', '\n epsilon')
                                        .replace('_mut_r', '\n mut_r')
                                        .replace('addcarbon', 'AddCarbon')
                                      for g in generators]
                        ax.set_yticklabels(new_labels, rotation=0, fontsize=14)
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
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    if save_name:
        if save_folder:
            plt.savefig(f'{save_folder}/{save_name}.svg', format="svg", bbox_inches='tight')
            plt.savefig(f'{save_folder}/{save_name}.png', format="png",dpi=300, bbox_inches='tight')
        else:
            plt.savefig(f'img/heat_map/{save_name}.svg', format="svg", bbox_inches='tight')
            plt.savefig(f'img/heat_map/{save_name}.png', format="png", dpi=300, bbox_inches='tight')
    plt.show()



def plot_combined_heatmap_with_single_column_for_each_metric_rotated(
        generators, receptors, scaffolds, splits,
        metrics=['TUPOR', 'SESY', 'ASER'], 
        title=None, save_name=None, using_norm_values=False,
        data_folder='', save_folder='',
        inter_metric_wspace=0.15,
        intra_metric_wspace=0.05
    ):

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

                    df = (preprocesing_org if not using_norm_values else preprocesing)(
                        type_cluster,
                        type_scaffold,
                        generators,
                        receptor,
                        data_folder
                    )

                    for met in metrics:
                        value = df[df.name.str.startswith(gen)][met].iloc[0]
                        data.append([
                            gen,
                            receptor,
                            type_scaffold,
                            type_cluster,
                            met,
                            value
                        ])

    df = pd.DataFrame(
        data,
        columns=['Generator', 'Receptor', 'Scaffold', 'Split', 'Metric', 'Value']
    )

    # -----------------------------------------
    # ROTATED LAYOUT
    # -----------------------------------------

    nrows = len(metrics)      # rows = metrics
    ncols = len(receptors)   # cols = receptors

    fig_width = 1.7 * 4 * ncols + 2
    fig_height = 6 * nrows * 1.3

    fig = plt.figure(figsize=(fig_width, fig_height))

    outer_gs = fig.add_gridspec(
        nrows=nrows,
        ncols=ncols,
        wspace=inter_metric_wspace,
        hspace=0.1
    )

    # MAIN LOOP
    for met_idx, metric in enumerate(metrics):

        metric_df = df[df['Metric'] == metric].copy()
        cmap_custom = make_cmap_to_white(metric_base_colors[metric])

        for rec_idx, receptor in enumerate(receptors):

            inner = outer_gs[met_idx, rec_idx].subgridspec(
                nrows=1,
                ncols=4,
                wspace=intra_metric_wspace,
                hspace=0.0
            )

            group_axes = []

            for sc_idx, scaffold_type in enumerate(["csk", "murcko"]):

                block_df = metric_df[
                    (metric_df['Receptor'] == receptor) &
                    (metric_df['Scaffold'] == scaffold_type)
                ]

                for split_idx, split in enumerate(["dis", "sim"]):

                    col = sc_idx * 2 + split_idx
                    ax = fig.add_subplot(inner[0, col])
                    group_axes.append(ax)

                    sub_df = block_df[block_df['Split'] == split].copy()
                    sub_df = sub_df.set_index('Generator').reindex(generators)

                    heatmap_array = sub_df['Value'].to_numpy().reshape(-1, 1)

                    vmin = heatmap_array.min()
                    vmax = heatmap_array.max()

                    heatmap_flat = heatmap_array.flatten()
                    max_idx = np.argmax(heatmap_flat)

                    annot_array = []
                    for i, val in enumerate(heatmap_flat):
                        txt = f"{val:.4f}" if not using_norm_values else f"{val:.3f}"
                        if i == max_idx:
                            txt = r"$\bf{" + txt + "}$"
                        annot_array.append(txt)

                    annot_array = np.array(annot_array).reshape(heatmap_array.shape)

                    show_colorbar = (sc_idx == 1 and split_idx == 1)

                    if show_colorbar:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes(
                            "right",
                            size="5%",
                            pad=0.05
                        )
                    else:
                        cax = None

                    sns.heatmap(
                        heatmap_array,
                        annot=annot_array,
                        fmt="",
                        cmap=cmap_custom,
                        ax=ax,
                        cbar=show_colorbar,
                        cbar_ax=cax,
                        annot_kws={
                            "size": 16,
                            "color": "black"
                        },
                        vmin=vmin,
                        vmax=vmax
                    )

                    ax.set_aspect("auto")

                    # X-axis = split labels
                    ax.set_xticks([0.5])
                    ax.set_xticklabels([split], rotation=0, ha="center", fontsize=14)

                    # -----------------------------------------
                    # Y LABELS  (GENS only in first column)
                    # -----------------------------------------
                    if rec_idx == 0 and sc_idx == 0 and split_idx == 0:
                        new_labels = [
                            g.replace('_epsilon', '\n epsilon')
                             .replace('_mut_r', '\n mut_r')
                             .replace('addcarbon', 'AddCarbon')
                            for g in generators
                        ]

                        ax.set_yticks(np.arange(len(generators)) + 0.5)
                        ax.set_yticklabels(new_labels, rotation=0, fontsize=15)
                        ax.set_ylabel(metric, fontsize=16, fontweight="bold", labelpad=8)

                    else:
                        ax.set_yticks([])
                        ax.set_ylabel("")

            # -----------------------------------------
            # TOP TITLE = RECEPTOR
            # -----------------------------------------
            if met_idx == 0:

                p0 = group_axes[0].get_position()
                p1 = group_axes[-1].get_position()

                x_mid = (p0.x0 + p1.x1) / 2
                y_top = p0.y1 + 0.018

                fig.text(
                    x_mid,
                    y_top,
                    receptor.replace('_', ' '),
                    ha="center",
                    va="bottom",
                    fontsize=16,
                    fontweight="bold"
                )

                # CSK = axes 0,1
                p0 = group_axes[0].get_position()
                p1 = group_axes[1].get_position()

                x_mid = (p0.x0 + p1.x1) / 2
                y_top = p0.y1 + 0.004

                fig.text(
                    x_mid,
                    y_top,
                    "CSK",
                    ha="center",
                    va="bottom",
                    fontsize=14
                )

                # MURCKO = axes 2,3
                p2 = group_axes[2].get_position()
                p3 = group_axes[3].get_position()

                x_mid = (p2.x0 + p3.x1) / 2
                y_top = p2.y1 + 0.004

                fig.text(
                    x_mid,
                    y_top,
                    "MURCKO",
                    ha="center",
                    va="bottom",
                    fontsize=14
                )

    # -----------------------------------------
    # FINAL TOUCHES
    # -----------------------------------------
    if title:
        fig.suptitle(title, fontsize=14, y=0.995)

    plt.tight_layout(rect=[0, 0, 1, 0.99])

    if save_name:
        if save_folder:
            plt.savefig(f'{save_folder}/{save_name}.svg', format="svg", bbox_inches='tight')
            plt.savefig(f'{save_folder}/{save_name}.png', format="png", dpi=300, bbox_inches='tight')
            plt.savefig(f'{save_folder}/{save_name}.pdf', bbox_inches='tight') 
        else:
            plt.savefig(f'img/heat_map/{save_name}.svg', format="svg", bbox_inches='tight')
            plt.savefig(f'img/heat_map/{save_name}.png', format="png", dpi=300, bbox_inches='tight')

    plt.show()
