import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
st.title('RMG input User Interface')
cif_file = st.file_uploader("upload cif file", type=".cif")
if cif_file:
  if not os.path.isdir("tempDir"):
    os.mkdir("tempDir")
  with open(os.path.join("tempDir", cif_file.name), "wb") as f:
    f.write(cif_file.getbuffer())
  filename = "tempDir/"+cif_file.name

  f = open(filename)
  st.write(f.read())
  rmginput = subprocess.check_output(["python","cif2rmg.py", filename])
  rmginput_str = rmginput.decode()
  #lines = rmginput_str.split("\n")
  #for line in lines:
  #  st.write(line)
  rmgfilename = cif_file.name.split(".")[0] +".rmg"
  st.download_button(
     label="Downlowd rmg input file",
     data= rmginput_str,
     file_name = rmgfilename)

