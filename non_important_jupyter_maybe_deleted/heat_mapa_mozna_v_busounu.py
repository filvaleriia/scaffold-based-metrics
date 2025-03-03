def plot_heatmaps_from_baseline(base_line, data_dict, type_split, scaf, receptor='', name_save=''):
    """
    Plots heatmaps comparing different subsets to a baseline.

    Args:
    - base_line (pd.DataFrame): DataFrame containing the baseline data to compare against.
    - data_dict (dict): Dictionary with keys representing subset names (e.g., '_10k', '_100k', '_500k', '') 
                        and values representing corresponding DataFrames.
    """

    num_subsets = len(data_dict)  # Dynamicky získáme počet subsetů
    fig, axes = plt.subplots(1, num_subsets, figsize=(8 * num_subsets, 8))  # 1 řádek, dynamický počet sloupců

    if num_subsets == 1:
        axes = [axes]  # Pokud je pouze jeden subset, přetvoříme axes na seznam pro kompatibilitu s iterací

    # Baseline pro porovnání
    baseline_df = base_line[['TUPOR', 'SESY', 'ASER']].copy()
    baseline_df.index = base_line.name.tolist()

    # Iterace přes subsety
    for k, (subset, df) in enumerate(data_dict.items()):
        # Načtení dat pro aktuální subset
        normalized_df = df[['TUPOR', 'SESY', 'ASER']]
        normalized_df.index = df.name.tolist()
        # Zajištění, že baseline a aktuální subset mají stejný index a sloupce
        baseline_df = baseline_df.set_index(normalized_df.index)
        ax = axes[k]  # Vybere subplot

        # Pokud je baseline, zobrazí všechny hodnoty; pro ostatní aplikuje masku
        if subset == '':
            mask = np.zeros_like(normalized_df, dtype=bool)  # Žádné maskování (zobrazí vše)
        else:
            mask = np.abs(normalized_df - baseline_df) < 0.1  # Maska pro hodnoty blízké baseline

        # Vykreslení heatmapy s maskou
        sns.heatmap(normalized_df, annot=True, cmap='viridis', cbar_kws={'label': 'Normalized Value'}, ax=ax, mask=mask)
        
        # Nastavení titulu pro subset
        ax.set_title(f'{scaf} - {subset if subset else "BASELINE"}', fontsize=16, wrap=True)

        # Změna názvů generátorů pro lepší čitelnost
        new_labels = [label.get_text().replace('_epsilon', '\n epsilon').replace('_mut_r', '\n mut_r') for label in ax.get_yticklabels()]
        ax.set_yticklabels(new_labels, rotation=0, ha="right", fontsize=14)

        # Nastavení popisků na ose X
        ax.set_xticks([0, 1, 2])  
        ax.set_xticklabels(['TUPOR', 'SESY', 'ASER'])

        ax.set_xticklabels(ax.get_xticklabels(), ha="center", fontsize=14)

    # Přidání celkového titulku
    fig.suptitle(f'Heatmaps with Differences from Baseline > 0.1 ({scaf}, {type_split}) Original values for {receptor}', fontsize=16)

    # Automatické přizpůsobení rozvržení
    plt.tight_layout(rect=[0, 0, 1, 0.95])  
    plt.savefig(f'img/heat_mapa/{receptor}/{name_save}.svg', format="svg")
    plt.show()





for type_cluster in ['dis', 'sim']:  # Iterate through different split types (dis, sim)
    for type_scaffold in ['csk', 'murcko']:  # Iterate through scaffold types (e.g., csk)
        receptor = 'Leukocyte_elastase'
        
        # Prepare a dictionary to store the data for different subsets
        # These subsets are differentiated by the suffixes (e.g., '_10k', '_100k', etc.)
        for subset in ['_10k', '_100k', '_500k', '']:
            generators_name_list = [
                f"Molpher{subset}",
                f"REINVENT{subset}",
                f"DrugEx_GT_epsilon_0.1{subset}",
                f"DrugEx_GT_epsilon_0.6{subset}",
                f"DrugEx_RNN_epsilon_0.1{subset}",
                f"DrugEx_RNN_epsilon_0.6{subset}",
                f"GB_GA_mut_r_0.01{subset}",
                f"GB_GA_mut_r_0.5{subset}",
                f"addcarbon{subset}"
            ]
            
            df = hv.preprocesing(type_cluster, type_scaffold, generators_name_list, receptor)
            
            # Load the baseline data for the given split type and scaffold
            if subset == '': #baseline
                base_line = df
                
            data_dict[subset] = df
        
        name_save = f'5x1_heat_maps_{type_cluster}_{type_scaffold}_{receptor}_diff_real'
        # Call the function 'plot_heatmaps_from_baseline' to plot the heatmaps for each subset and compare them with the baseline
        hv.plot_heatmaps_from_baseline(base_line, data_dict, type_cluster, type_scaffold, receptor = receptor, name_save = name_save)
