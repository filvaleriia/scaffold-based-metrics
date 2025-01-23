import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
    page_title="Metrics",
    page_icon="📝",
)
st.sidebar.header("Metrics")


st.title('Metrics')

st.markdown(
    """
    #### **TUPOR**
    > this is the main metric. It reflects the main goal of chemical structural generators, which is to generate new biological active compounds. So, this metric is about recovery of active scaffolds. 

    #### **ASER** 
    > this metric reflects biological activity too. This metric is about how many compounds generator found with active scaffolds. 

    #### **SESY** 
    > this metric is about how diverse our virtual library is. It is important for showing which generator has more different compounds in the output.

    | NAME | FORMULA | LIMITATIONS |
    |-----------------|-|-|
    |▪ True positive recall all (TUPOR)| $$TUPOR = {{UASo} \over UASr}$$ | [0 ; 1] |
    |▪ Set scaffold yield (SESY)| $$SESY = {{USo} \over SSo}$$ | [0 ; 1] |
    |▪ Absolute set scaffold recall (ASER)| $$ASER = {{CwASo} \over SSo}$$ |[0 ; 1] |
    |▪ Active scaffold ratio (ASR)| $$ASR = {{UASo} \over {{USo} - {UASo}}}$$ | [0 ; 1] |
 

    * UASo - unique active scaffolds in cOS
    * USo - unique scaffolds in cOS (all unique scaffold = unique active scaffolds + unique other scaffolds)
    * UASr - unique active scaffolds in cRS
    * USr - unique scaffolds in cRS
    * SSo - cOS size 
    * CwASo - Compounds with Active Scaffol in cOS 


 
"""
)