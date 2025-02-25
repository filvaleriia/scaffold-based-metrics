import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def normalize_and_plot_radar_charts(df, type_split, type_scaffold, subset, receptor):
    """
    This function processes data for different scaffold types and splits, 
    normalizes the values, and creates radar charts to compare different generators.
    The radar charts are saved as SVG files and displayed.
    """

    # Loop through different data splits and scaffold types
    
    # Define the labels for the metrics to be normalized
    labels = ['TUPOR', 'SESY', 'ASER', 'ACR']
    num_metrics = len(labels)
    
    # Normalize values using min-max scaling (0-1 range)
    ranges = {}  # Dictionary to store min and max values for each label
    for label in labels:
        min_val = df[label].min()
        max_val = df[label].max()
        df[label] = (df[label] - min_val) / (max_val - min_val)  # Normalize the column
        ranges[label] = (min_val, max_val)  # Store the range for each label
    
    # Define the generator names, labels, and styles for plotting
    generators = [
        {'name': f'Molpher{subset}_mean', 'label': 'Molpher', 'linestyle': '-', 'marker': 'o', 'color': 'blue'},
        {'name': f'REINVENT{subset}_mean', 'label': 'REINVENT', 'linestyle': '--', 'marker': 's', 'color': 'green'},
        {'name': f'DrugEx_GT_epsilon_0.1{subset}_mean', 'label': 'DrugEx_GT_epsilon_0.1', 'linestyle': '-.', 'marker': 'D', 'color': 'purple'},
        {'name': f'DrugEx_GT_epsilon_0.6{subset}_mean', 'label': 'DrugEx_GT_epsilon_0.6', 'linestyle': ':', 'marker': 'v', 'color': 'orange'},
        {'name': f'DrugEx_RNN_epsilon_0.1{subset}_mean', 'label': 'DrugEx_RNN_epsilon_0.1', 'linestyle': '-.', 'marker': '*', 'color': 'magenta'},
        {'name': f'DrugEx_RNN_epsilon_0.6{subset}_mean', 'label': 'DrugEx_RNN_epsilon_0.6', 'linestyle': '--', 'marker': '<', 'color': 'brown'},
        {'name': f'GB_GA_mut_r_0.01{subset}_mean', 'label': 'GB_GA_mut_r_0.01', 'linestyle': ':', 'marker': '>', 'color': 'gray'},
        {'name': f'GB_GA_mut_r_0.5{subset}_mean', 'label': 'GB_GA_mut_r_0.5', 'linestyle': ':', 'marker': '>', 'color': 'black'},
        {'name': f'addcarbon{subset}_mean', 'label': 'AddCarbon', 'linestyle': '-', 'marker': 'h', 'color': 'pink'},
    ]
    
    # Create angles for the radar chart axes
    angles = np.linspace(0, 2 * np.pi, num_metrics, endpoint=False).tolist()
    angles += [angles[0]]  # Close the circle
    
    # Create the figure and axis for plotting
    fig, ax = plt.subplots(figsize=(13, 10), subplot_kw=dict(polar=True))
    
    # Plot each generator
    for gen in generators:
        generator_name = gen['name']
        label = gen['label']
        linestyle = gen['linestyle']
        marker = gen['marker']
        color = gen['color']
        
        # Get the values for the generator
        values = df[df.name == generator_name][labels].values.flatten().tolist()
        if not values:  # Skip if no data is found for the generator
            print(f"No data found for {generator_name}")
            continue
        values += [values[0]]  # Close the circle
        
        # Plot the generator's data on the radar chart
        ax.plot(angles, values, label=label, linestyle=linestyle, marker=marker, color=color)
        ax.fill(angles, values, alpha=0.25, color=color)
    
    # Adjust the radar chart settings
    ax.set_theta_offset(np.pi / 2)  # Set the top of the chart to be the starting point
    ax.set_theta_direction(-1)  # Reverse the direction of the chart
    
    # Set the axis labels with ranges
    xticks_labels = [f"{label}\n({ranges[label][0]:.2f}-{ranges[label][1]:.2f})" for label in labels]
    #xticks_labels = [f"{label}" for label in labels]
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(xticks_labels)
    ax.tick_params(pad=20)
    
    # Set the value range for the axes (0-1, as values are normalized)
    ax.set_ylim(0, 1)
    
    # Add the legend
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.05))
    
    # Add the title
    if subset == '':
        subset_text = 'base'
    else:
        subset_text = subset
    plt.title(f'Comparison of Generators Based on Normalized Metrics {type_scaffold} {type_split} {subset_text}', size=16)
    
    # Save the plot as an SVG file
    plt.savefig(f'img/radar_plots/{receptor}/radar_results_all_{type_scaffold}_{type_split}{subset}.svg', format="svg")
    
    # Show the plot
    plt.show()

