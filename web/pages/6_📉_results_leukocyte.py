import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
    page_title="Results for Leukocyte elastase",
    page_icon="📉",
)
st.sidebar.header("Results for Leukocyte elastase:")

st.sidebar.markdown('**Molpher**')
st.sidebar.markdown('**DrugEx GT**')
st.sidebar.markdown('**DrugEx RNN**')
st.sidebar.markdown('**REINVENT**')
st.sidebar.markdown('**GB_GA**')
st.sidebar.markdown('**Add carbone**')
st.sidebar.markdown('**VAE generator**')

st.subheader('Molpher:')
st.markdown(
    """
    Molpher generator with Dissimilarity split and CSK scaffolds
"""
)

generator = 'Molpher'
receptor = 'Leukocyte_elastase'
type_scaffold = 'csk'
type_split = 'dis'

df = pd.read_csv(f"../data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
st.dataframe(data = df)

st.markdown(
    """
    Molpher generator with Similarity split and CSK scaffolds
"""
)

if st.checkbox('Show statistics about runs for dissimilarity split'):
    st.markdown(
    """
    Statistic about good runs for informations, but because I made one mistakes with input data I did some Ad hoc changes. So this statistic it's only for showing how many runs is ended for GLucocorticoid receptor
    """
    )
    df = pd.read_csv(f"../data/information_about_runs/{receptor}/{type_split}/information_about_runs.csv")
    st.dataframe(data = df)




# SIMILARITY SPLIT
type_split = 'sim'

df = pd.read_csv(f"../data/results/{receptor}/{type_scaffold}_scaffolds/{type_split}/{generator}/df_all_clusters_with_mean.csv")
st.dataframe(data = df)

if st.checkbox('Show statistics about runs for similarity split'):
    st.markdown(
    """
    Statistic about good runs for informations, but because I made one mistakes with input data I did some Ad hoc changes. So this statistic it's only for showing how many runs is ended for GLucocorticoid receptor
    """
    )
    df = pd.read_csv(f"../data/information_about_runs/{receptor}/{type_split}/information_about_runs.csv")
    st.dataframe(data = df)
