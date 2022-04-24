
import streamlit as st
def add_pseudo(species):  
  expand_ = st.expander("pseudopotential")
  with expand_:
      st.subheader('Pseudopotentials')
      cstart, col1 = st.columns([0.1,1])
      pseudo_dict={}
      with col1:
          pp_oncv = st.checkbox("ONCV: build-in normal conserving Haman pseudopotentials", True)
          pp_gbrv = st.checkbox("GBRV: build-in ultra-soft pseudopotentials", False) 
            
          pseudolines = ""  
          if pp_oncv:
              pseudolines = 'internal_pseudo_type = "nc_accuracy"\n'
          if pp_gbrv:
              pseudolines = 'internal_pseudo_type = "ultrasoft"\n'
          select_pp = st.checkbox("select pseudopotentials files by yourself", False)

      cstart, col1 = st.columns([0.2,1])
      with col1:
          if select_pp:
             pseudo_dir = st.text_input("pseudopotential file directory", value="./", on_change=None)
             pseudolines = 'pseudo_dir = "' + pseudo_dir + '"\n'
             pseudolines += 'pseudopotential = "\n' 
             for sp in species:
                 pseudo_dict[sp] = st.text_input(sp+":", value=sp+".UPF", on_change=None)
                 pseudolines += sp + '   ' + pseudo_dict[sp] +'\n'

             pseudolines += '"\n' 
      return pseudolines       
