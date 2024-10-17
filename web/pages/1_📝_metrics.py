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

    |             | |
    |-------------|-|
    |▪ True positive recall all (TUPOR) – | $$TUPOR = {{UASo} \over UASr}$$ |
    |▪ True positive recall all (TUPORunique) – | $$TUPORunique = {{UASo} \over (USo * UASr)}$$ |
    |▪ Set scaffold yield (SESY) –| $$SESY = {{USo} \over SSo}$$ |
    |▪ Absolute set scaffold recall (ASER) –| $$ASER = {{CwASo} \over SSo}$$ |


    * UASo - unique active scaffolds in cOS
    * USo - unique scaffolds in cOS (all unique scaffold = unique active scaffolds + unique other scaffolds)
    * UASr - unique active scaffolds in cRS
    * SSo - cOS size 
    * CwASo - Compounds with Active Scaffol in cOS 


 
"""
)