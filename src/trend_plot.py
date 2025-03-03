import matplotlib.pyplot as plt
import pandas as pd

def preprocesing(type_cluster, type_scaffold, generators_name_list, receptor):
    # Define path to data
    link = f"data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}"
    link_mean = [f"{link}/{generator}/{generator}_mean_{type_scaffold}_{type_cluster}.csv" for generator in generators_name_list]
    
    # Load data
    df_list = [pd.read_csv(f) for f in link_mean]
    df = pd.concat(df_list, axis=0, ignore_index=True)

    df["name"] = df["name"].str.replace("_mean", "", regex=False)

    
    return df

def plot_trends(df, type_split, type_scaf, subset, receptor, name_save):
    """
    Generates two trend plots from the given DataFrame:
    - One for 'TUPOR' and 'SESY'
    - One for 'ASER' an
    
    Parameters:
    - df (pd.DataFrame): Data containing 'name', 'TUPOR', 'SESY', 'ASER', an columns.
    - type_split (str): The type of data split (e.g., 'dis' or 'sim').
    - type_scaf (str): The scaffold type (e.g., 'csk' or 'murcko').
    - subset (str): The subset identifier (e.g., '_10k', '_100k', '_500k', or '').
    """
    df = df[['name', 'TUPOR', 'SESY', 'ASER']]
    # Ensure required columns exist in DataFrame
    required_columns = {'name', 'TUPOR', 'SESY', 'ASER'}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"Missing required columns in DataFrame. Expected columns: {required_columns}")

    # Create two subplots side by side
    fig, axes = plt.subplots(1, 1, figsize=(16, 6))  # 1 row, 2 columns
    
    # Left plot for TUPOR and SESY
    axes.plot(df["name"], df["TUPOR"], marker='o', label="TUPOR")
    axes.plot(df["name"], df["SESY"], marker='o', label="SESY")
    axes.plot(df["name"], df["ASER"], marker='o', label="ASER")
    axes.set_title(f"TUPOR and SESY trends ({type_split}, {type_scaf})", fontsize=14)
    axes.set_xlabel("Generator", fontsize=12)
    axes.set_ylabel("Value of metrics", fontsize=12)

    x_labels = df['name'].str.replace('_mean', '').tolist()
    x_positions = range(len(x_labels)) 

    axes.set_xticks(x_positions)  
    axes.set_xticklabels(x_labels, rotation=45, ha='right')
    axes.legend(title="Metrics")

    
    # Adjust layout
    plt.subplots_adjust(bottom=0.25)
    plt.tight_layout()
    fig.suptitle(f'Trends plots for {type_split} split and {type_scaf} scaffold for {subset} subset', y=1.02, fontsize = 16) 
    plt.savefig(f'img/trends_plots/{receptor}/{name_save}.svg', format="svg")
    # Show the plot
    plt.show()



def plot_subset_trends(type_scaf, df_all, subsets_list, receptor, name_save):
    """
    Generates two trend plots for a given scaffold type:
    - One for 'TUPOR' and 'SESY'
    - One for 'ASER' an
    
    Parameters:
    - type_scaf (str): The scaffold type ('csk' or 'murcko').
    - subsets_list (list of str): The list of subsets to use (e.g., ['_10k', '_100k', '_500k', '']).
    """
    
    # Initialize an empty DataFrame to store combined data
   

    # Create the plots
    fig, axes = plt.subplots(2, 1, figsize=(18, 12))  # Two subplots in one column
    plt.style.use("seaborn-darkgrid")
    
    # Define metrics for each plot
    metrics_1 = ['TUPOR', 'SESY']
    metrics_2 = ['ASER']
    
    # Define colors and styles
    subset_styles = ['-', ':']  # Line styles for subsets
    colors = ['b', 'g', 'r', 'c', 'm']  # Colors for different splits
    
    # First plot: TUPOR and SESY
    for i, subset in enumerate(subsets_list):
        subset_label = subset if subset else 'base'
        for k, type_split in enumerate(['dis', 'sim']):
            for j, col in enumerate(metrics_1):
                df_subset = df_all[(df_all['Subset'] == subset_label) & (df_all['Split'] == type_split)]
                axes[0].plot(
                    df_subset['name'],
                    df_subset[col],
                    label=f"{col} {type_split} ({subset_label})",
                    linestyle=subset_styles[j % len(subset_styles)],
                    color=colors[k % len(colors)],
                    marker='o'
                )

    axes[0].set_title(f"Trends for TUPOR and SESY ({type_scaf})", fontsize=16)
    axes[0].set_ylabel("Value of metrics", fontsize=14)
    axes[0].legend(title="Metrics and Subsets", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)

    x_labels = df_all['name'].unique()
    x_positions = range(len(x_labels)) 

    axes[0].set_xticks(x_positions)  
    axes[0].set_xticklabels(x_labels, rotation=60, ha='right')

    # Second plot: ASER and ACR
    colors = ['r', 'black']
    for i, subset in enumerate(subsets_list):
        subset_label = subset if subset else 'base'
        for k, type_split in enumerate(['dis', 'sim']):
            for j, col in enumerate(metrics_2):
                df_subset = df_all[(df_all['Subset'] == subset_label) & (df_all['Split'] == type_split)]
                axes[1].plot(
                    df_subset['name'].str.replace('_mean', '').tolist(),
                    df_subset[col],
                    label=f"{col} {type_split} ({subset_label})",
                    linestyle=subset_styles[j % len(subset_styles)],
                    color=colors[k % len(colors)],
                    marker='o'
                )

    axes[1].set_title(f"Trends for ASER and ACR ({type_scaf})", fontsize=16)
    axes[1].set_ylabel("Value of metrics", fontsize=14)
    axes[1].legend(title="Metrics and Subsets", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    axes[1].set_xticks(x_positions)  
    axes[1].set_xticklabels(x_labels, rotation=60, ha='right')

    # Adjust layout
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f'img/trends_plots/{receptor}/{name_save}.svg', format="svg")
    # Show the plots
    plt.show()


def plot_trends_base(scaffolds=['csk', 'murcko'], type_cluster=['dis', 'sim'], df_all = pd.DataFrame(), receptor='', name_save=''):
    """
    Function to plot trends for given scaffolds and type_cluster.
    
    Arguments:
    scaffolds -- List of scaffold types to consider (default is ['csk', 'murcko'])
    type_cluster -- List of type_cluster to consider (default is ['dis', 'sim'])
    """
    
    # Set up the figure and axes for plotting
    fig, axes = plt.subplots(2, 2, figsize=(20, 15))
    plt.style.use("seaborn-darkgrid")
    
    # Metrics for each graph
    metrics = [['TUPOR', 'SESY'], ['ASER']]
    colors = ['b', 'r', 'g', 'black']
    
    # Loop over the metrics for plotting
    for k, metrics_1 in enumerate(metrics):
        row = k
        if k == 0:
            subset_styles = ['-', ':']
        else:
            subset_styles = ['-', ':']
        
        # Loop over scaffold types
        for i, scaffold in enumerate(scaffolds):
            col = i % 2
            l = 0
            # Loop over split types (DIS, SIM)
            for k, split in enumerate(type_cluster):
                subset_style_num = 0
                # Loop over metrics for plotting
                for col_metric in metrics_1:
                    df_subset = df_all[(df_all['Scaffold'] == scaffold) & (df_all['Type_Split'] == split)]
                    generators = df_subset['name'].str.replace('_mean', '').tolist()
                    
                    # Add error bars (if needed)
                    # error_bars = [std_dict[f'{gen}_{scaffold}_{split}'][metrics_1.index(col_metric)] for gen in generators]
                    
                    # Plot the data with error bars
                    axes[row, col].errorbar(
                        generators, df_subset[col_metric],
                        # yerr=error_bars,  # Add error bars
                        label=f"{col_metric} ({scaffold} {split})",
                        linestyle=subset_styles[subset_style_num],
                        color=colors[l % len(colors)],
                        marker='o',
                        # capsize=5  # Set the size of the error bar caps
                    )
                    subset_style_num += 1
                l += 1
            axes[row, col].set_title(f"Trends for {' and '.join(metrics_1)} ({scaffold})", fontsize=12)
            axes[row, col].set_ylabel("Value of metrics", fontsize=10)
            axes[row, col].legend(title="Metrics and Scaffold + Split", fontsize=8)
            
            # Set the positions for x-axis labels
            axes[row, col].set_xticks(range(len(generators)))  # Set positions for X
            axes[row, col].set_xticklabels(generators, rotation=60, fontsize=10, ha='right')
    
    # Adjust layout to avoid overlap
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f'img/trends_plots/{receptor}/{name_save}.svg', format="svg")
    plt.show()