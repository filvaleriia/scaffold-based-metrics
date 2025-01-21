import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
    page_title="Results for Glucocorticoid receptor",
    page_icon="📈",
)
st.sidebar.header("Results for Glucocorticoid receptor:")

st.sidebar.markdown('**Molpher**')
st.sidebar.markdown('**DrugEx GT**')
st.sidebar.markdown('**DrugEx RNN**')
st.sidebar.markdown('**REINVENT**')
st.sidebar.markdown('**GB_GA**')
st.sidebar.markdown('**Add carbone**')
st.sidebar.markdown('**VAE generator**')

pd.options.display.float_format = '{:12.5e}'.format

web = True
if web == True:
    current_directory = '/mount/src'
else:
    current_directory = '/workspaces'
st.header("Show results for Glucocorticoid receptor:")
receptor = 'Glucocorticoid_receptor'

st.markdown(f'##### Results for Glucocorticoid receptor  dissimilartity split')
if st.checkbox(f'Show mean results for dissimilarity split: csk scaffold first then murcko to compare results between different scaf. First table then img and then next results'):
    type_split = 'dis'
    
    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
            df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
            st.dataframe(data = df)


            col1, col2 = st.columns(2)
            image_width = 300
            # V prvním sloupci zobrazíme první obrázek
            with col1:
                st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

            # Ve druhém sloupci zobrazíme druhý obrázek
            with col2:
                st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)


if st.checkbox(f'Show mean results for dissimilarity split: csk scaffold first then murcko to compare results between different scaf. First all tables then all imgs '):
    type_split = 'dis'

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
            df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
            st.dataframe(data = df)

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            col1, col2 = st.columns(2)
            image_width = 300
            # V prvním sloupci zobrazíme první obrázek
            with col1:
                st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

            # Ve druhém sloupci zobrazíme druhý obrázek
            with col2:
                st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)

#SIMILARITY SPLIT
st.markdown(f'##### Results for Glucocorticoid receptor  similartity split')
if st.checkbox(f'Show mean results for similarity split: csk scaffold first then murcko to compare results between different scaf. First table then img and then next results'):
    type_split = 'sim'
    

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
            df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
            st.dataframe(data = df)


            col1, col2 = st.columns(2)
            image_width = 300
            # V prvním sloupci zobrazíme první obrázek
            with col1:
                st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

            # Ve druhém sloupci zobrazíme druhý obrázek
            with col2:
                st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)


if st.checkbox(f'Show mean results for similarity split: csk scaffold first then murcko to compare results between different scaf. First all tables then all imgs '):
    type_split = 'sim'
    

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
            df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
            st.dataframe(data = df)

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            col1, col2 = st.columns(2)
            image_width = 300
            # V prvním sloupci zobrazíme první obrázek
            with col1:
                st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

            # Ve druhém sloupci zobrazíme druhý obrázek
            with col2:
                st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)

#RESULTS TO COMPARE CSK/MURCKO DISSIMILARITY SPLIT
st.markdown(f'##### Results to compare csk and murcko scaffolds')

if st.checkbox(f'Show mean results for compare csk/murcko scaf for dissimilatiry split: csk scaffold first then murcko to compare results between different scaf. First table then img and then next results '):
    type_split = 'dis'
    for _ in ['', '_1k', '_10k', '_100k', '_500k']:
        for type_split in ['dis']:
            for type_scaffold in ['csk', 'murcko']:
            
        
                #compare_differen_results----------------------------------------------------

                df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
                df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
                st.dataframe(data = df)
            
            for type_scaffold in ['csk', 'murcko']:
                col1, col2 = st.columns(2)
                image_width = 300
                # V prvním sloupci zobrazíme první obrázek
                with col1:
                    st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

                # Ve druhém sloupci zobrazíme druhý obrázek
                with col2:
                    st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)

#RESULTS TO COMPARE CSK/MURCKO SIMILARITY SPLIT
if st.checkbox(f'Show mean results for compare csk/murcko scaf for similatiry split: csk scaffold first then murcko to compare results between different scaf. First table then img and then next results '):
    for _ in ['', '_1k', '_10k', '_100k', '_500k']:
        for type_split in ['sim']:
            for type_scaffold in ['csk', 'murcko']:
            
        
                #compare_differen_results----------------------------------------------------

                df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
                df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
                st.dataframe(data = df)
            
            for type_scaffold in ['csk', 'murcko']:
                col1, col2 = st.columns(2)
                image_width = 300
                # V prvním sloupci zobrazíme první obrázek
                with col1:
                    st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

                # Ve druhém sloupci zobrazíme druhý obrázek
                with col2:
                    st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)



#RESULTS to compare dis/sim
st.markdown(f'##### Results to compare dissimilarity split and similarity split')
if st.checkbox(f'Show mean results for compare dis/sim: csk scaffold first then murcko to compare results between different scaf. First table then img and then next results '):
    type_split = 'dis'
    for _ in ['', '_1k', '_10k', '_100k', '_500k']:
        for type_scaffold in ['csk', 'murcko']:
            for type_split in ['dis', 'sim']:
        
                #compare_differen_results----------------------------------------------------

                df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
                df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
                st.dataframe(data = df)
            
            for type_split in ['dis', 'sim']:
                col1, col2 = st.columns(2)
                image_width = 300
                # V prvním sloupci zobrazíme první obrázek
                with col1:
                    st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

                # Ve druhém sloupci zobrazíme druhý obrázek
                with col2:
                    st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)


st.markdown(f'##### Results for Glucocorticoid receptor  dissimilartity split without subsets')
if st.checkbox(f'Show mean results for dissimilarity split without subsets: csk scaffold first then murcko to compare results between different scaf. First table then img and then next results'):
    type_split = 'dis'
    
    for type_scaffold in ['csk', 'murcko']:
        for _ in ['']:
            #compare_differen_results----------------------------------------------------

            df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
            df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
            st.dataframe(data = df)


            col1, col2 = st.columns(2)
            image_width = 300
            # V prvním sloupci zobrazíme první obrázek
            with col1:
                st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

            # Ve druhém sloupci zobrazíme druhý obrázek
            with col2:
                st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)


if st.checkbox(f'Show mean results for dissimilarity split without subsets: csk scaffold first then murcko to compare results between different scaf. First all tables then all imgs '):
    type_split = 'dis'

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
            df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
            st.dataframe(data = df)

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            col1, col2 = st.columns(2)
            image_width = 300
            # V prvním sloupci zobrazíme první obrázek
            with col1:
                st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

            # Ve druhém sloupci zobrazíme druhý obrázek
            with col2:
                st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)

#SIMILARITY SPLIT
st.markdown(f'##### Results for Glucocorticoid receptor  similartity split without subsets')
if st.checkbox(f'Show mean results for similarity split without subsets: csk scaffold first then murcko to compare results between different scaf. First table then img and then next results'):
    type_split = 'sim'
    

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
            df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
            st.dataframe(data = df)


            col1, col2 = st.columns(2)
            image_width = 300
            # V prvním sloupci zobrazíme první obrázek
            with col1:
                st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

            # Ve druhém sloupci zobrazíme druhý obrázek
            with col2:
                st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)


if st.checkbox(f'Show mean results for similarity split without subsets: csk scaffold first then murcko to compare results between different scaf. First all tables then all imgs '):
    type_split = 'sim'
    

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}{_}.csv")
            df = df[['name','type_cluster','scaffold', 'TUPOR', 'SESY', 'ASER','ASR']]
            st.dataframe(data = df)

    for type_scaffold in ['csk', 'murcko']:
        for _ in ['', '_1k', '_10k', '_100k', '_500k']:
            #compare_differen_results----------------------------------------------------

            col1, col2 = st.columns(2)
            image_width = 300
            # V prvním sloupci zobrazíme první obrázek
            with col1:
                st.image(f'{current_directory}/recall_metrics/img/results/new_plots/radar_results_all_{type_scaffold}_{type_split}_{_}.png', width=image_width, use_column_width=True)

            # Ve druhém sloupci zobrazíme druhý obrázek
            with col2:
                st.image(f'{current_directory}/recall_metrics/img/heat_mapa/heat_mapa_{type_split}_{type_scaffold}_{_}.png', width=image_width, use_column_width=True)

st.markdown(f'##### Results for Glucocorticoid receptor for individual generators')
if st.checkbox(f'Show mean results for individual generators each generator individualy'):
    
    for generator in ['Molpher', 'REINVENT', 'DrugEx_GT_epsilon_0.1', 'DrugEx_GT_epsilon_0.6', 'DrugEx_RNN_epsilon_0.1', 'DrugEx_RNN_epsilon_0.6', 'GB_GA_mut_r_0.01', 'addcarbon']:
        if st.checkbox(f'{generator}'):
        
            st.markdown(f"#### {generator} ")
            for type_split in ['dis', 'sim']:
                for type_scaffold in ['csk', 'murcko']:
                
                    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
                    df = df[['name','type_cluster','scaffold','TUPOR_', 'TUPOR', 'SESY', 'ASER', 'ASR']]
                    st.dataframe(data = df)

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
                    
                        st.markdown(f"#### {generator} ")


                        df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
                        df = df[['name','type_cluster','scaffold','TUPOR_', 'TUPOR', 'SESY', 'ASER', 'ASR']]
                        st.dataframe(data = df)

                        if type_scaffold == 'csk' and generator == 'Molpher': #showing only once
                        
                            if st.checkbox(f'Show statistics about runs for dissimilarity split {type_split}'):
                                st.markdown(
                                """
                                Statistic about good runs for informations, but because I made one mistakes with input data I did some Ad hoc changes. So this statistic it's only for showing how many runs is ended for GLucocorticoid receptor
                                """
                                )
                                df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_runs/{receptor}/{type_split}/information_about_runs.csv")
                                st.dataframe(data = df)
   


st.markdown(f"## t-SNE vizualization for each generator")
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
