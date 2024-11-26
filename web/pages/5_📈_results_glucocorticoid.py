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

    #DrugEx_GT_1-----------------------------------------------------------
    generator = 'DrugEx_GT_1'

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
    #pustit vsehcny metriky na similarity!

    st.markdown(f"#### Compare mean of metrics value between different generators")

    df = pd.read_csv(f"{current_directory}/recall_metrics/data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/mean_{type_scaffold}_{type_split}.csv")
    #df = df[['name','scaffold', 'USo', 'SSo', 'TUPOR', 'TUPOR_2', 'SESY', 'ASER', 'CwASo']]
    st.dataframe(data = df)

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

    #DrugEx_GT_1-----------------------------------------------------------
    generator = 'DrugEx_GT_1'

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

    

st.markdown(f"##t-SNE vizualization for each generator")
st.markdown(
    """
    t-SNE for dis 0
"""
)
st.image(f'{current_directory}/recall_metrics/hist_all_dis_0.png', caption='Glucocorticoid_receptor')
