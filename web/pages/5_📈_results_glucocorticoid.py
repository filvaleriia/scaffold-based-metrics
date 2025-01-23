import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
    page_title="Results for Glucocorticoid receptor",
    page_icon="📈",
)
st.sidebar.header("Results for Glucocorticoid Receptor:")
generators = [
    "Molpher", "DrugEx GT", "DrugEx RNN", 
    "REINVENT", "GB_GA", "Add carbone", "VAE generator"
]
for generator in generators:
    st.sidebar.markdown(f"**{generator}**")

pd.options.display.float_format = '{:12.5e}'.format

web = True
if web == True:
    current_directory = '/mount/src'
else:
    current_directory = '/workspaces'

#POMOCNE FUNKCE
def display_data_and_images(type_split, scaffold_type, subset, display_images=True, display_table = True):
    file_path = f"{current_directory}/recall_metrics/data/results/{receptor}/{scaffold_type}_scaffolds/{type_split}/mean_{scaffold_type}_{type_split}{subset}.csv"
    img1_path = f"{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{scaffold_type}_{type_split}_{subset}.png"
    img2_path = f"{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{scaffold_type}_{subset}.png"
    
    # Zobrazeni tabulek
    if display_table:
        df = pd.read_csv(file_path)
        df = df[["name", "type_cluster", "scaffold", "TUPOR", "SESY", "ASER", "ASR"]]
        st.dataframe(df)

    
    # Zobrazení obrázků (pokud je potřeba)
    if display_images:
        col1, col2 = st.columns(2)
        with col1:
            st.image(img1_path)
        with col2:
            st.image(img2_path)


#Hlavi telo
st.header("Show results for Glucocorticoid receptor:")
receptor = 'Glucocorticoid_receptor'

# Sekce: Dissimilarity splitu
st.markdown("### Dissimilarity Split Results")
if st.checkbox(f'Show mean results for dissimilarity split: (Table-img results)'):
    type_split = 'dis'  
    for type_scaffold in ['csk', 'murcko']:
        for subset in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=True)  # Tabulka - obrazek


if st.checkbox(f'Show mean results for dissimilarity split: (All tables then all img) '):
    type_split = 'dis'
    for type_scaffold in ['csk', 'murcko']:
        for subset in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky

    for type_scaffold in ['csk', 'murcko']:
        for subset in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=False)  # Obrazky

#SEKCE: SIMILARITY SPLIT
st.markdown("### Similarity Split Results")
if st.checkbox(f'Show mean results for similarity split: (Table-img results)'):
    type_split = 'sim'
    for type_scaffold in ['csk', 'murcko']:
        for subset in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=True)  # Tabulka - obrazek


if st.checkbox(f'Show mean results for similarity split: (All tables then all img)'):
    type_split = 'sim'
    for type_scaffold in ['csk', 'murcko']:
        for subset in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky

    for type_scaffold in ['csk', 'murcko']:
        for subset in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=False)  # Obrazky

#RESULTS TO COMPARE CSK/MURCKO DISSIMILARITY SPLIT
st.markdown("### Comparison of CSK and Murcko for Dissimilarity Split")

if st.checkbox("Compare CSK/Murcko scaffolds (dissimilarity split):"):
    type_split = 'dis'
    for subset in ['', '_1k', '_10k', '_100k', '_500k']:
        for type_split in ['dis']:
            for type_scaffold in ['csk', 'murcko']:
            
                display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky
            
            for type_scaffold in ['csk', 'murcko']:
                display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=False)  # Obrazky

#RESULTS TO COMPARE CSK/MURCKO SIMILARITY SPLIT
st.markdown("### Comparison of CSK and Murcko for Similarity Split")

if st.checkbox("Compare CSK/Murcko scaffolds (similarity split):"):
    for subset in ['', '_1k', '_10k', '_100k', '_500k']:
        for type_split in ['sim']:
            for type_scaffold in ['csk', 'murcko']:
                display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky
            
            for type_scaffold in ['csk', 'murcko']:
                display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=False)  # Obrazky



#RESULTS to compare dis/sim
st.markdown("### Comparison of Dissimilarity and Similarity Splits")
if st.checkbox("Compare dissimilarity and similarity splits"):
    type_split = 'dis'
    for subset in ['', '_1k', '_10k', '_100k', '_500k']:
        for type_scaffold in ['csk', 'murcko']:
            for type_split in ['dis', 'sim']:
        
                display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky
            
            for type_split in ['dis', 'sim']:
                display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=False)  # Obrazky


st.markdown(f'### Comparison of CSK and Murcko for Dissimilarity Split without subsets')
if st.checkbox(f'Show mean results for dissimilarity split without subsets: (Table-img results)'):
    type_split = 'dis'
    for type_scaffold in ['csk', 'murcko']:
        for subset in ['']:
            display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=True)  # Tabulka-obrazka


if st.checkbox(f'Show mean results for dissimilarity split without subsets: (All tables then all img)'):
    type_split = 'dis'

    for type_scaffold in ['csk', 'murcko']:
        for subset in ['']:
            display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky

    for type_scaffold in ['csk', 'murcko']:
        for subset in ['']:
            display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=False)  # Obrazky

#SIMILARITY SPLIT
st.markdown(f'### Comparison of CSK and Murcko for Similarity Split without subsets')
if st.checkbox(f'Show mean results for similarity split without subsets: (Table-img results)'):
    type_split = 'sim'
    for type_scaffold in ['csk', 'murcko']:
        for subset in ['']:
            display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=True)  # Tabulka-obrazka


if st.checkbox(f'Show mean results for similarity split without subsets: (All tables then all img)'):
    type_split = 'sim'

    for type_scaffold in ['csk', 'murcko']:
        for subset in ['']:
            display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky

    for type_scaffold in ['csk', 'murcko']:
        for subset in ['']:
            display_data_and_images(type_split, type_scaffold, subset, display_images=True, display_table=False)  # Obrazky

st.markdown("### Results for Individual Generators")
if st.checkbox("Show results for each generator individually"):
    
    for generator in ['Molpher', 'REINVENT', 'DrugEx_GT_epsilon_0.1', 'DrugEx_GT_epsilon_0.6', 'DrugEx_RNN_epsilon_0.1', 'DrugEx_RNN_epsilon_0.6', 'GB_GA_mut_r_0.01', 'addcarbon']:
        if st.checkbox(f'{generator}'):
            st.markdown(f"#### {generator} ")
            for type_split in ['dis', 'sim']:
                for type_scaffold in ['csk', 'murcko']:
                
                    display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky

                    if type_scaffold == 'csk' and generator == 'Molpher': #showing only once
                    
                        if st.checkbox(f'Show statistics about runs for dissimilarity split {type_split}'):
                            st.markdown(
                            """
                            Statistic about good runs for informations, but because I made one mistakes with input data I did some Ad hoc changes. So this statistic it's only for showing how many runs is ended for GLucocorticoid receptor
                            """
                            )
                            df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_runs/{receptor}/{type_split}/information_about_runs.csv")
                            st.dataframe(data = df)


if st.checkbox(f'Show mean results for individual to compare between generators'):
    for type_split in ['dis', 'sim']:
        if st.checkbox(f'{type_split}'):
            for type_scaffold in ['csk', 'murcko']:
                if st.checkbox(f'{type_scaffold}'):
                    for generator in ['Molpher', 'REINVENT', 'DrugEx_GT_epsilon_0.1', 'DrugEx_GT_epsilon_0.6', 'DrugEx_RNN_epsilon_0.1', 'DrugEx_RNN_epsilon_0.6', 'GB_GA_mut_r_0.01', 'addcarbon']:
                        subset = ''
                        st.markdown(f"#### {generator} ")
                        display_data_and_images(type_split, type_scaffold, subset, display_images=False, display_table=True)  # Tabulky

                        if type_scaffold == 'csk' and generator == 'Molpher': #showing only once
                        
                            if st.checkbox(f'Show statistics about runs for dissimilarity split {type_split}'):
                                st.markdown(
                                """
                                Statistic about good runs for informations, but because I made one mistakes with input data I did some Ad hoc changes. So this statistic it's only for showing how many runs is ended for GLucocorticoid receptor
                                """
                                )
                                df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_runs/{receptor}/{type_split}/information_about_runs.csv")
                                st.dataframe(data = df)
   

#ZOBRAZOVANI OVERLAPU:
st.markdown(f"## Overlaps plots ")
if st.checkbox(f'Overplaps plots all'):
    st.image(f"{current_directory}/recall_metrics/img/overlaps/triples/combined_overlaps_csk.png")


#ZOBRAZOVANI OVERLAPU:
st.markdown(f"## Box plots ")

for type_scaf in ['csk', 'murcko']:
    if st.checkbox(f'Box plots for {type_scaf}'):
        for type_split in ['dis','sim']:
            st.image(f"{current_directory}/recall_metrics/img/box_plots/{type_scaf}/combined_box_plot_{type_split}.png")


#ZOBRAZOVANI TRENDU:
st.markdown(f"## Trends plot ")
if st.checkbox(f'Trends_plot'):
    if st.checkbox(f'Show subsets comparison'):
        for type_scaf in ['csk', 'murcko']:
            for type_split in ['dis','sim']:
                st.image(f"{current_directory}/recall_metrics/img/trands/{type_scaf}/trends_combined_all_subsets_{type_split}.png")
                st.image(f"{current_directory}/recall_metrics/img/trands/{type_scaf}/trends_combined_adjusted_subsets_{type_split}.png")
    if st.checkbox(f'Show split/scaf comparison'):
        st.image(f"{current_directory}/recall_metrics/img/trands/trends_combined_TUPOR_SESY.png")
        st.image(f"{current_directory}/recall_metrics/img/trands/trends_combined_ASER_ASR.png")

#ZOBRAZOVANI Heat mapa:
st.markdown(f"## Heat map ")
if st.checkbox(f'Heat maps 1x5'):
    
    for type_scaf in ['csk', 'murcko']:
        for type_split in ['dis','sim']:
            st.image(f"{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaf}_1x5.png")

if st.checkbox(f'Heat map for base subset'):
    st.image(f"{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_base_all.png")




#ZOBRAZOVANI T-SNE
st.markdown(f"## t-SNE vizualization for each generator")
if st.checkbox(f't-SNE vizualization for each generator'):
    for number in [0,1,2,3,4]:
        if st.checkbox(f'Show t-SNE vizualization for dissimilarity split {number}'):
            type_cluster = 'dis'

            for generator in ['addcarbon', 'Molpher', 'REINVENT', 'GB_GA', 'DrugEx_GT',
                              'DrugEx_GT_epsilon_0.1', 'DrugEx_GT_epsilon_0.6', 'DrugEx_RNN',
                              'DrugEx_RNN_epsilon_0.1', 'DrugEx_RNN_epsilon_0.6']:
                st.image(f'{current_directory}/recall_metrics/img/t-SNE/{generator}/t_sne_gen_is_rs_{type_cluster}_{number}.png')

    for number in [0,1,2,3,4]:
        if st.checkbox(f'Show t-SNE vizualization for similarity split {number}'):
            type_cluster = 'sim'

            for generator in ['addcarbon', 'Molpher', 'REINVENT', 'GB_GA', 'DrugEx_GT',
                              'DrugEx_GT_epsilon_0.1', 'DrugEx_GT_epsilon_0.6', 'DrugEx_RNN',
                              'DrugEx_RNN_epsilon_0.1', 'DrugEx_RNN_epsilon_0.6']:
                st.image(f'{current_directory}/recall_metrics/img/t-SNE/{generator}/t_sne_gen_is_rs_{type_cluster}_{number}.png')

    st.markdown(f"## t-SNE vizualization for DrugEx with different epsilon")
    for number in [0,1,2,3,4]:
        if st.checkbox(f'Show t-SNE vizualization DrugEx for dissimilarity split {number}'):
            type_cluster = 'dis'

            for generator in ['DrugEx_epsilon_0.01', 'DrugEx_epsilon_0.1', 'DrugEx_epsilon_0.2',
                                'DrugEx_epsilon_0.4', 'DrugEx_epsilon_0.6']:
                st.image(f'{current_directory}/recall_metrics/img/t-SNE/DrugEx/{generator}_{type_cluster}_{number}.png')

    for number in [0,1,2,3,4]:
        if st.checkbox(f'Show t-SNE vizualization DrugEx for similarity split {number}'):
            type_cluster = 'sim'

            for generator in ['DrugEx_epsilon_0.01', 'DrugEx_epsilon_0.1', 'DrugEx_epsilon_0.2',
                                'DrugEx_epsilon_0.4', 'DrugEx_epsilon_0.6']:
                st.image(f'{current_directory}/recall_metrics/img/t-SNE/DrugEx/{generator}_{type_cluster}_{number}.png')

st.markdown(f"## Percentage statistics of compounds with desirable properties")
if st.checkbox(f'Show statistics'):
    st.image(f'{current_directory}/recall_metrics/img/sa_qsprpred/QSPRPRED_SA_1.png', caption='Molpher,REINVENT,addcarbon')
    st.image(f'{current_directory}/recall_metrics/img/sa_qsprpred/QSPRPRED_SA_2.png', caption='DrugEx_GT')
    st.image(f'{current_directory}/recall_metrics/img/sa_qsprpred/QSPRPRED_SA_3.png', caption='DrugEx_RNN')
