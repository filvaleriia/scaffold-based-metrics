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
    |▪ True positive recall all (TUPOR) – | $$TUPOR = {{NAS} \over UAS}$$ |
    |▪ Set scaffold yield (SESY) –| $$SESY = {{NS} \over SS}$$ |
    |▪ Absolute set scaffold recall (ASER) –| $$ASER = {{tRS} \over SS}$$ |


    * tRS - count of active compounds which have active scaffold in cOS
    * SS - cOS size 
    * NS - unique new scaffolds in cOS
    * NAS - unique active scaffolds in cRS
    * UAS - unique scaffods in cRS 

 
"""
)