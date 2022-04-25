import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
from rmg_parser import *
from add_items import *
st.title('RMG input User Interface')
st.write('<style>div.row-widget.stRadio > div{flex-direction:row;justify-content: left}<style>',
        unsafe_allow_html=True)

uploaded_file = st.file_uploader("Uploda a file")
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
  #st.write(crmg.species)
  description = st.text_input("description", value="description of the input file")
  rmginput_str = 'description="'+description+'"  \n'

  grid_lines = add_grid(crmg.cell)
  pseudo_lines = add_pseudo(crmg.species)
  kpoint_lines = add_kpoints(crmg.cell)
  ctrl_lines = add_control()
  scf_lines = add_scf()
  mixing_lines = add_mixing()
  xc_lines = add_xc()

      
  rmginput_str += grid_lines
  rmginput_str += scf_lines
  rmginput_str += mixing_lines
  rmginput_str += xc_lines
  rmginput_str += ctrl_lines
  rmginput_str += kpoint_lines
  rmginput_str += pseudo_lines
  rmginput_str += crmg.rmginput
  rmgfilename = os.path.basename(filename).split(".")[0] +".rmg"
  st.download_button(
     label="Downlowd rmg input file",
     data=rmginput_str,
     file_name = rmgfilename)
  show_rmginput = st.checkbox("show the generated rmg input file", False)
  if show_rmginput:
    st.markdown(rmginput_str)

