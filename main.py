import streamlit as st
import pandas as pd
import numpy as np
import subprocess
st.title('RMG input User Interface')

rmginput = subprocess.check_output(["python","a.py"])
st.write( rmginput)

