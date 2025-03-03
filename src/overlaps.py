import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import pandas as pd
import random
from itertools import combinations
from pathlib import Path
import os

# Function to draw the Venn diagram with percentages
def draw_venn_multiple(sets_dict, subset_keys, ax):
    subset_dict = {k: sets_dict[k] for k in subset_keys}
    subsets = [subset_dict[key] for key in subset_keys]
    total_unique = len(set.union(*subsets))

    # If there is no data to plot (empty sets), don't plot and just hide the axis
    if total_unique == 0:
        ax.set_visible(False)
        return

    # Calculate the sizes of overlaps
    counts = {
        '100': len(subsets[0] - subsets[1] - subsets[2]),
        '010': len(subsets[1] - subsets[0] - subsets[2]),
        '001': len(subsets[2] - subsets[0] - subsets[1]),
        '110': len(subsets[0] & subsets[1] - subsets[2]),
        '101': len(subsets[0] & subsets[2] - subsets[1]),
        '011': len(subsets[1] & subsets[2] - subsets[0]),
        '111': len(subsets[0] & subsets[1] & subsets[2]),
    }

    # Convert the counts to percentages
    counts_percentage = {k: (v / total_unique * 100 if total_unique > 0 else 0) for k, v in counts.items()}

    # Draw the Venn diagram
    venn = venn3(subsets=[subsets[0], subsets[1], subsets[2]], set_labels=subset_keys, alpha=0.5, ax=ax)

    # Set percentage values on the diagram
    for region, value in counts_percentage.items():
        label = venn.get_label_by_id(region)
        if label:
            label.set_text(f"{value:.1f}%")

    # Customize generator labels
    for label in venn.set_labels:
        label.set_fontsize(10)
        label.set_text(label.get_text().replace("_epsilon", "\n epsilon"))

# Function to process the data and plot the Venn diagrams
def plot_overlaps(receptor, generators, num_graphs):
    scaffold_sets = {}
    main_dir = Path(__file__).resolve().parents[1]  # Main directory for data and images

    # Load the scaffold sets
    for type_scaf in ['csk', 'murcko']:
        # Load data into sets
        for generator in generators:
            scaffold_sets[generator] = set()
            for num in range(5):  # Cluster numbers 0-4
                file_path = f"{main_dir}/data/results/{receptor}/{type_scaf}_scaffolds/dis/{generator}/scaffolds_of_output_set_cluster_{num}_dis_{generator}.csv"
                try:
                    df = pd.read_csv(file_path, header=None)
                    scaffold_sets[generator].update(df[0].tolist())
                except FileNotFoundError:
                    print(f"File {file_path} not found!")

        # Select the specified number of random generator combinations
        all_combinations = list(combinations(generators, 3))
        if num_graphs > len(all_combinations):
            print(f"Requested number of graphs is greater than the available combinations. Selecting all {len(all_combinations)} combinations.")
            selected_combinations = all_combinations
        else:
            selected_combinations = random.sample(all_combinations, num_graphs)

        # If there are no combinations to display, skip plotting
        if not selected_combinations:
            print(f"No combinations available to display. Skipping plotting.")
            continue

        # Dynamically calculate the number of rows
        rows = (num_graphs // 4) + (num_graphs % 4 > 0)  # Adjust rows based on graphs count
        fig, axes = plt.subplots(rows, 4, figsize=(20, 5 * rows))  # Adjust height dynamically
        axes = axes.flatten()

        # Plot the selected overlap graphs
        for i, comb in enumerate(selected_combinations):
            draw_venn_multiple(scaffold_sets, comb, axes[i])

        # Hide unused axes
        for j in range(i+1, len(axes)):
            axes[j].set_visible(False)

        # Adjust layout
        plt.tight_layout()
        fig.suptitle(f'Overlaps Visualizations for {receptor} ({type_scaf})', fontsize=20)
        # Save the figure with a name based on the scaffold type
        output_folder = f'{main_dir}/img/overlaps/{receptor}'
        os.makedirs(output_folder, exist_ok=True)  # Ensure the folder exists
        plt.savefig(f'{output_folder}/selected_overlaps_{type_scaf}_{num_graphs}_graphs.png', dpi=500)

        # Display the graphs
        plt.show()

# Example usage:
# plot_overlaps('Leukocyte_elastase', ['Molpher', 'REINVENT', 'DrugEx_GT_epsilon_0.1', 'DrugEx_RNN_epsilon_0.1'], 12)
