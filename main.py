import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
from rmg_parser import *
from add_pseudo import *
st.title('RMG input User Interface')
st.write('<style>div.row-widget.stRadio > div{flex-direction:row;justify-content: center}<style>',
        unsafe_allow_html=True)
uploaded_file = st.file_uploader("Upload a file")
example_ =  st.checkbox("use an example cif file FeAs.cif", False)

filetype_set = [".cif", ".xyz"]
if uploaded_file:
  if not os.path.isdir("tempDir"):
    os.mkdir("tempDir")
  with open(os.path.join("tempDir", uploaded_file.name), "wb") as f:
    f.write(uploaded_file.getbuffer())
  filename = "tempDir/"+uploaded_file.name
  name_split = os.path.splitext(filename)
  if len(name_split) >1: filext = name_split[len(name_split)-1]
  if filext in filetype_set: 
      filetype = filext
  else:
      filetype = st.text_input("filetype")
elif example_:
  filename = "cifs/FeAs.cif"  
  filetype = ".cif"

if uploaded_file or example_:
  crmg = rmg_interface(filename, filetype)
  rmginput_str = crmg.rmginput
  #st.write(crmg.species)
  pseudolines = add_pseudo(crmg.species)

      



  rmginput_str += pseudolines
  rmgfilename = os.path.basename(filename).split(".")[0] +".rmg"
  st.download_button(
     label="Downlowd rmg input file",
     data=rmginput_str,
     file_name = rmgfilename)
  show_rmginput = st.checkbox("show the generated rmg input file", False)
  if show_rmginput:
    st.markdown(rmginput_str)

