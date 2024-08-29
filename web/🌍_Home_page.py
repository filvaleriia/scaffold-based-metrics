import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
    page_title="Recall Metrics project",
    page_icon="🌍",
)
st.sidebar.header("Recall Metrics project")


st.title('Recall Metrics project')

st.markdown(
    """

    The main goal of this project is to propose better benchmarks for comparing various generators. The motivation is that a sensible benchmark should capture properties of compounds that are interesting from a biological standpoint and not just "trivial" properties such as validity or uniqueness. The main such property, for which generators are deployed, is biological activity. The problem with biological activity is that QSAR models have their applicability domain, which complicates their use on entirely new chemotypes, whose design using molecular generators is, on the other hand, very desirable. For comparing generators, pharmacophore modeling or molecular docking could be used, which do not suffer from the problem of the applicability domain. However, in our project, we go a different way. Our metrics are based on scaffold recall, i.e., we obtain a set of compounds and their scaffolds which we know are biologically active, part of this set is used as input for the generator (its "compound Input Set") and in the generated virtual library, we look for scaffolds from the Recall Set, i.e., scaffolds that belong to active substances but were not used for training the generator. Thus, our metrics capture the ability of the generator to generate new active chemotypes. 

    In this project, we proposed three metrics to evaluate the performance of chemical structure generators and have the options to compare one generator with another. 

    We have chosen one (or two) targets to demonstrate the functionality of these metrics. 

    ### Important informations: 


    - ChEMBL 31: download 18.1.2023
    - conda environment is rdkit-env(conda activate rdkit-env)
    - data are download in 22.02.2023
    - Used descriptors: Morgan Fingerprint r=3, 2048
    - Threshold for nuclear receptor - <0-100> nM
    - Threshold for protease - <0-1000> nM
    - IC50 - pouzita merena aktivita
    - for splitting to 5 clusters used Kmedoids algorithm with Tanimoto distance (because in K-medoids we can change the distance, because in K-means they used Euclidian distance)

    ### Important shortcuts

    - cIS = compound Input Set
    - sIS = scaffold Input Set
    - cRS = compound Recall Set
    - sRS = scaffold Recall Set
    - cOS = compound Output Set
    - dis = Dissimilarity split
    - sim = Similarity split

"""
)