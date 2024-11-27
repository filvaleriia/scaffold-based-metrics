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

web = False
if web == True:
    current_directory = '/mount/src'
else:
    current_directory = '/workspaces'

for type_scaffold in ['csk', 'murcko']:

    st.header(f'Results for Glucocorticoid receptor used {type_scaffold} scaffold')

    # DISSIMILARITY SPLIT---------------------------------------------------------------------------------------------
    receptor = 'Glucocorticoid_receptor'
    #type_scaffold = 'csk'
    type_split = 'dis'
    st.markdown('### * Dissimilarity split')

    #compare_differen_results-----------------------------------------------------------


    st.markdown(f"#### Compare mean of metrics value between different generators")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}.csv")
    df = df[['name','scaffold', 'USo', 'SSo', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    if st.checkbox(f'Show vizualization of results dissimilarity split {type_scaffold} scaffold'):
                
        st.image(f'{current_directory}/recall_metrics/img/results/radar_results_{type_scaffold}_{type_split}.png')


    #MOLPHER--------------------------------------------------------
    generator = 'Molpher'

    st.markdown(f"#### {generator} ")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    if type_scaffold == 'csk': #showing only once

        if st.checkbox('Show statistics about runs for dissimilarity split'):
            st.markdown(
            """
            Statistic about good runs for informations, but because I made one mistakes with input data I did some Ad hoc changes. So this statistic it's only for showing how many runs is ended for GLucocorticoid receptor
            """
            )
            df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_runs/{receptor}/{type_split}/information_about_runs.csv")
            st.dataframe(data = df)

    #DrugEx_GT-----------------------------------------------------------
    generator = 'DrugEx_GT'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_GT_epsilon_0.1-----------------------------------------------------------
    generator = 'DrugEx_GT_epsilon_0.1'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_GT_epsilon_0.6-----------------------------------------------------------
    generator = 'DrugEx_GT_epsilon_0.6'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_RNN-----------------------------------------------------------
    generator = 'DrugEx_RNN'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_RNN_epsilon_0.1-----------------------------------------------------------
    generator = 'DrugEx_RNN_epsilon_0.1'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_RNN_epsilon_0.6-----------------------------------------------------------
    generator = 'DrugEx_RNN_epsilon_0.6'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)



    #GB_GA-----------------------------------------------------------
    generator = 'GB_GA'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)


    #add_carbon-----------------------------------------------------------
    generator = 'addcarbon'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)






    # SIMILARITY SPLIT---------------------------------------------------------------------------------------------------

    st.markdown('### * Similarity split')
    type_split = 'sim'

    #compare_differen_results-----------------------------------------------------------

    st.markdown(f"#### Compare mean of metrics value between different generators")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}.csv")
    df = df[['name','scaffold', 'USo', 'SSo', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    if st.checkbox(f'Show vizualization of results similarity split {type_scaffold} scaffold'):
                
        st.image(f'{current_directory}/recall_metrics/img/results/radar_results_{type_scaffold}_{type_split}.png')


    #MOLPHER
    generator = 'Molpher'

    st.markdown(f"#### {generator} ")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    if type_scaffold == 'csk': #showing only once
        if st.checkbox('Show statistics about runs for similarity split'):
            st.markdown(
            """
            Statistic about good runs for informations, but because I made one mistakes with input data I did some Ad hoc changes. So this statistic it's only for showing how many runs is ended for GLucocorticoid receptor
            """
            )
            df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_runs/{receptor}/{type_split}/information_about_runs.csv")
            st.dataframe(data = df)

    #DrugEx_GT-----------------------------------------------------------
    generator = 'DrugEx_GT'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_GT_epsilon_0.1-----------------------------------------------------------
    generator = 'DrugEx_GT_epsilon_0.1'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_GT_epsilon_0.6-----------------------------------------------------------
    generator = 'DrugEx_GT_epsilon_0.6'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_RNN-----------------------------------------------------------
    generator = 'DrugEx_RNN'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_RNN_epsilon_0.1-----------------------------------------------------------
    generator = 'DrugEx_RNN_epsilon_0.1'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

    #DrugEx_RNN_epsilon_0.6-----------------------------------------------------------
    generator = 'DrugEx_RNN_epsilon_0.6'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)


    #GB_GA-----------------------------------------------------------
    generator = 'GB_GA'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)


    #add_carbon-----------------------------------------------------------
    generator = 'addcarbon'

    st.markdown(f"#### {generator}")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
    df = df[['name','scaffold', 'USo', 'SSo','TUPOR_', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
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