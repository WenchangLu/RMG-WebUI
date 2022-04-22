import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
from cif2rmg import *
st.title('RMG input User Interface')
cif_file = st.file_uploader("upload cif file", type=".cif")
if cif_file:
  if not os.path.isdir("tempDir"):
    os.mkdir("tempDir")
  with open(os.path.join("tempDir", cif_file.name), "wb") as f:
    f.write(cif_file.getbuffer())
  filename = "tempDir/"+cif_file.name

  crmg = cifrmg_interface()
  rmginput_str = crmg.cif2rmg_run(filename)

  rmgfilename = cif_file.name.split(".")[0] +".rmg"
  st.download_button(
     label="Downlowd rmg input file",
     data=rmginput_str,
     file_name = rmgfilename)
  show_rmginput = st.checkbox("show the generated rmg input file", False)
  if show_rmginput:
    lines = rmginput_str.split("\n")
    for line in lines:
      st.write(line)
  #st.markdown(rmginput_str)
  #st.text(rmginput_str)
