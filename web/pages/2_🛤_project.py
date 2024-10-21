import streamlit as st
import pandas as pd
import numpy as np
from PIL import Image
import os

st.set_page_config(
    page_title="Project",
    page_icon="🛤",
)
st.sidebar.header("Project")

current_directory = os.getcwd()

# Zobrazení aktuálního adresáře ve Streamlit aplikaci
st.write(f"Current working directory: {current_directory}")


st.title('Project')

st.markdown(
    """
    The goal of the project is to experimentally test the behavior of metrics I have designed for two different molecular generators and compare them with each other. Within the scope of this project, two molecular generators will be tested: DrugEx, which is based on deep learning, and Molpher, which is based on evolutionary algorithms.
"""
)

uploaded_file = st.file_uploader(
    "../img/main_steps.jpg", accept_multiple_files=False)
if uploaded_file is not None:
    file_name = uploaded_file
else:
    file_name = "DatabaseSample.xlsx"

st.image('/mount/src/recall_metrics/img/main_steps.jpg', caption='Main steps of this project')

st.subheader('Split into Input Set and Recall Set:')
st.markdown(
    """
    For splitting data to clusters we convert data to CSK scaffolds and try to split data based the scaffold similarity(tanimoto distance) using K-medoids algorithm
    We divided all data for 5 clusters, and then for "dissimilarity" part we selectd four clusters for sIS and another one to sRS. And we done this five times. Than for creating cIS and cRS we found all compound wich have scaffold in sIS and sRS.
"""
)

st.image('/workspaces/recall_metrics/img/dissimilarity_split.jpg', caption='Dissimilarity split')

st.markdown(
    """
We divided all data for 5 clusters, and then for "similarity" part we select 20% from each clusters for sRS and other 80% from each cluster for sIS.
"""
)

st.image('/workspaces/recall_metrics/img/similarity_split.jpg', caption='Similarity split')

