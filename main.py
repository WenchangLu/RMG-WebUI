import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
from cif2rmg import *
st.title('RMG input User Interface')
st.write('<style>div.row-widget.stRadio > div{flex-direction:row;justify-content: center}<style>',
        unsafe_allow_html=True)
cif_file = st.file_uploader("Upload a CIF file",type=".cif")
example_ =  st.checkbox("use an example cif file FeAs.cif", False)
if cif_file:
  if not os.path.isdir("tempDir"):
    os.mkdir("tempDir")
  with open(os.path.join("tempDir", cif_file.name), "wb") as f:
    f.write(cif_file.getbuffer())
  filename = "tempDir/"+cif_file.name
elif example_:
  filename = "cifs/FeAs.cif"  

if cif_file or example_:
  crmg = cifrmg_interface()
  rmginput_str = crmg.cif2rmg_run(filename)
  #st.write(crmg.species)
  st.subheader('Pseudopotentials')
  pseudo_dict={}
  pp_oncv = st.checkbox("ONCV: build-in normal conserving Haman pseudopotentials", True)
  pp_gbrv = st.checkbox("GBRV: build-in ultra-soft pseudopotentials", False) 
    
  pseudolines = ""  
  if pp_oncv:
      pseudolines = 'internal_pseudo_type = "nc_accuracy"\n'
  if pp_gbrv:
      pseudolines = 'internal_pseudo_type = "ultrasoft"\n'

  select_pp = st.checkbox("select pseudopotentials files by yourself", False)
  if select_pp:
      pseudo_dir = st.text_input("pseudopotential file directory", value="./", on_change=None)
      pseudolines = 'pseudo_dir = "' + pseudo_dir + '"\n'
      pseudolines += 'pseudopotential = "\n' 
      for sp in crmg.species:
          pseudo_dict[sp] = st.text_input(sp+":", value=sp+".UPF", on_change=None)
          pseudolines += sp + '   ' + pseudo_dict[sp] +'\n'

      pseudolines += '"\n' 

      



  rmginput_str += pseudolines
  rmgfilename = os.path.basename(filename).split(".")[0] +".rmg"
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

