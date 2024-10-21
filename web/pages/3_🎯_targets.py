import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
    page_title="Target",
    page_icon="🎯",
)
st.sidebar.header("Targets:")

web = True
if web == True:
    current_directory = '/mount/src'
else:
    current_directory = '/workspaces'

st.sidebar.markdown('**Glucocorticoid receptor**')
st.sidebar.markdown('**Leukocyte elastase**')

st.subheader('Glucocorticoid receptor:')
receptor = 'Glucocorticoid_receptor'
st.markdown(
    """
    The glucocorticoid receptor (GR) is a type of nuclear receptor that binds glucocorticoids, which are a class of steroid hormones. This receptor plays a critical role in regulating a wide range of physiological processes, including metabolism, immune response, and stress response
"""
)
st.image(f'{current_directory}/recall_metrics/img/glucocorticoid_receptor.png', caption='Glucocorticoid_receptor')


if st.checkbox('Show unique scaffolds in each cluster for Glucocorticoid receptor'):
    st.subheader('Unique CSK scaffolds for each cluster')
    for x in range(5):
        st.image(f'{current_directory}/recall_metrics/img/clusters/scaffolds_in_cluster_{x}_Glucocorticoid_receptor.png', caption=f'CSK scaffolds for cluster {x}')

if st.checkbox('Show statistics about clusters for Glucocorticoid_receptor'):
    st.subheader('Mean value between clusters')
    st.markdown("This statistic created in the follow step: create all combinations between one csk in one cluster and another csk in the second clusters, then I calculated Morgan Fingerprints with R=10, nBits = 2048, beacuse this statistics based on CSK scaffolds, so for better representation bits")
    df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_clusters/{receptor}/mean_value_between_clusters.csv")
    st.dataframe(data = df)

    st.subheader('Min value between clusters')
    df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_clusters/{receptor}/min_value_between_clusters.csv")
    st.dataframe(data = df)

    st.subheader('Max value between clusters')
    df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_clusters/{receptor}/max_value_between_clusters.csv")
    st.dataframe(data = df)

    st.subheader('Statistic inside clusters: mean, min, max values')
    df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_clusters/{receptor}/statistic_incide_cluster.csv")
    st.dataframe(data = df)

#Leukocyte_elastase-----------------------------------------------------------------------------------------------
st.subheader('Leukocyte elastase:')
receptor = 'Leukocyte_elastase'
st.markdown(
    """
    Leukocyte elastase, also known as neutrophil elastase, is a serine protease enzyme primarily found in neutrophils, a type of white blood cell. It plays a crucial role in the body's immune response by degrading various proteins, thus helping to control infections and inflammation. However, its activity must be tightly regulated, as excessive or uncontrolled elastase activity can lead to tissue damage and contribute to various diseases.
"""
)

st.image(f'{current_directory}/recall_metrics/img/leukocyte_elastase.png', caption='Dissimilarity split')


if st.checkbox('Show unique scaffolds in each cluster for Leukocute elastase'):
    st.subheader('Unique CSK scaffolds for each cluster')
    for x in range(5):
        st.image(f'{current_directory}/recall_metrics/img/clusters/scaffolds_in_cluster_{x}_Leukocyte_elastase.png', caption=f'CSK scaffolds for cluster {x}')

if st.checkbox('Show statistics about clusters for Lukocyte elastase'):

    st.markdown("This statistic created in the follow step: create all combinations between one csk in one cluster and another csk in the second clusters, then I calculated Morgan Fingerprints with R=10, nBits = 2048, beacuse this statistics based on CSK scaffolds, so for better representation bits")

    st.subheader('Mean value between clusters')
    df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_clusters/{receptor}/mean_value_between_clusters.csv")
    st.dataframe(data = df)

    st.subheader('Min value between clusters')
    df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_clusters/{receptor}/min_value_between_clusters.csv")
    st.dataframe(data = df)

    st.subheader('Max value between clusters')
    df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_clusters/{receptor}/max_value_between_clusters.csv")
    st.dataframe(data = df)

    st.subheader('Statistic inside clusters: mean, min, max values')
    df = pd.read_csv(f"{current_directory}/recall_metrics/data/information_about_clusters/{receptor}/statistic_incide_cluster.csv")
    st.dataframe(data = df)