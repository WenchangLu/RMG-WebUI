
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
              pseudolines = 'internal_pseudo_type = "nc_accuracy"  \n'
          if pp_gbrv:
              pseudolines = 'internal_pseudo_type = "ultrasoft"  \n'
          select_pp = st.checkbox("select pseudopotentials files by yourself", False)

      cstart, col1 = st.columns([0.2,1])
      with col1:
          if select_pp:
             pseudo_dir = st.text_input("pseudopotential file directory", value="./", on_change=None)
             pseudolines = 'pseudo_dir = "' + pseudo_dir + '"  \n'
             pseudolines += 'pseudopotential = "  \n' 
             for sp in species:
                 pseudo_dict[sp] = st.text_input(sp+":", value=sp+".UPF", on_change=None)
                 pseudolines += sp + '   ' + pseudo_dict[sp] +'  \n'

             pseudolines += '"  \n' 
      return pseudolines    
def add_kpoint_mesh(cell):
    cs, col1, col2, col3 = st.columns([0.2,1,1,1])
    with col1:
        k_delta_str = st.text_input("kdelta(2PI/bohr)", value="0.2", help ="use kdelta to estimate kmesh")
    k_delta = float(k_delta_str)
    recip_lat= cell.reciprocal_latticevectors()
    for i in range(3):
        for j in range(3):
            recip_lat[i][j] = recip_lat[i][j]/cell.lengthscale
    kmesh_init = [max(1, int(b.length()/k_delta)) for b in recip_lat]

    kmesh_init_str =""
    for i in range(3):
        kmesh_init_str += str(kmesh_init[i]) + "  "
    with col2:
        kmesh_str = st.text_input("kpoint mesh", value=kmesh_init_str)
    with col3:
        kshift_str = st.text_input("kpoint shift", value="0 0 0", help="0 0 0 including Gamma point")
    kpointlines = 'kpoint_mesh="' + kmesh_str +'"  \n'    
    kpointlines += 'kpoint_shift="' + kshift_str +'"  \n'    
    return kpointlines

          

def add_kpoints(cell):
    expand_ = st.expander("choose k points")
    with expand_:
        kpointlines = add_kpoint_mesh(cell)
    return kpointlines
