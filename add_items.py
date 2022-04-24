
import streamlit as st
def add_pseudo(species):  
  expand_ = st.expander("PSEUDOPOTENTIAL")
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
    kmesh_init = [max(1, int(round(b.length()/k_delta))) for b in recip_lat]

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
    expand_ = st.expander("CHOOSE K POINTS")
    with expand_:
        kpointlines = add_kpoint_mesh(cell)
    return kpointlines
def add_control():
    expand_ = st.expander("CONTROL OPTIONS")
    with expand_:
        start_mode = st.radio("start mode", 
                ["LCAO Start", "Restart From File", "Random Start",
                  "FIREBALL Start", "Gaussian Start",
                   "Start TDDFT", "Restart TDDFT",
                    "Modified LCAO Start"])
        calculation_mode= st.radio("calculation mode", 
                ["Quench Electron  ", 
                    "Relax Structure  ", 
                    "Constant Volume And Energy  ",
                    "Constant Temperature And Energy   ",
                    "Constant Pressure And Energy  ", 
                    "Plot  ", 
                    "Psi Plot  ", 
                    "Band Structure Only  ", 
                    "NEB Relax  ",
                    "Dimer Relax  ",
                    "TDDFT "])

        kohn_sham_solver=st.radio("kohn_sham_solver", 
                ["davidson", "multigrid"],
                help="Davidson is prefered for a small system and multigrid for a large system")

        subdiag_driver = st.radio("diagonalizatoin libs",
                ["auto", "lapack", "scalapack", "magma", 
                 "cusolver", "elpa", "rocsolver"])

        vdw_corr = st.radio("van der Waals correction", 
                ["None", "DFT-D2", "Grimme-D2","DFT-D3"])
        vdwdf_grid = st.radio("grid for vdw corr",
                ["Coarse", "Fine"])
        relax_mass = st.radio("mass for atoms", ["Atomic", "Equal"], 
                help="equal mas for fast relax may help in some cases")
        dos_method = st.radio("density of state calc", 
              ["tetrahedra", "Gaussian"])

        occupations_type = st.radio("occupation type",
                ["Fermi Dirac", "Fixed", "Cold Smearing", "MethfesselPaxton"])
        occ_smear = st.text_input("occupation smear in eV", value ="0.1")

        poisson_solver = st.radio("Poisson Solver",
                ["pfft", "multigrid"])

        md_tem_ctrl = st.radio("MD temperature control", 
                ["Nose Hoover Chains","Anderson Rescaling"])
        md_integration_order = st.radio("MD Integration order",
                ["2nd Velocity Verlet",
                 "3rd Beeman-Velocity Verlet",
                 "5th Beeman-Velocity Verlet"])
        ctrl_lines = ""
        ctrl_lines += 'start_mode          ="' +start_mode +'"  \n'
        ctrl_lines += 'calculation_mode    ="' +calculation_mode +'"  \n'
        ctrl_lines += 'kohn_sham_solver    ="' +kohn_sham_solver +'"  \n'
        ctrl_lines += 'subdiag_driver      ="' +subdiag_driver +'"  \n'
        ctrl_lines += 'vdw_corr            ="' +vdw_corr +'"  \n'
        ctrl_lines += 'vdwdf_grid          ="' +vdwdf_grid +'"  \n'
        ctrl_lines += 'relax_mass          ="' +relax_mass +'"  \n'
        ctrl_lines += 'dos_method          ="' +dos_method +'"  \n'
        ctrl_lines += 'occupations_type    ="' +occupations_type +'"  \n'
        ctrl_lines += 'occupation_electron_temperature_eV="' +occ_smear +'"  \n'
        ctrl_lines += 'poisson_solver      ="' +poisson_solver +'"  \n'
        ctrl_lines += 'md_tem_ctrl         ="' +md_tem_ctrl +'"  \n'
        ctrl_lines += 'md_integration_order="' +md_integration_order +'"  \n'
        return ctrl_lines

def add_grid(cell):
    expand_ = st.expander("CHOOSE REAL SPACE GRID")
    with expand_:
        cs, col1, col2, col3 = st.columns([0.1,1,2,1])
        grid_spacing_str = col1.text_input("grid spacing(bohr)", value="0.35",
                    help ="use grid spacing to determine the real space grid")
        grid_spacing = float(grid_spacing_str)
        nx = int(round(cell.a/grid_spacing))
        ny = int(round(cell.b/grid_spacing))
        nz = int(round(cell.c/grid_spacing))
        i2 = 1
        for i in range(4):
            i2 *= 2
            nx1 = (nx+i2-1)//i2 * i2
            ny1 = (ny+i2-1)//i2 * i2
            nz1 = (nz+i2-1)//i2 * i2
            h_max = max(cell.a/nx1, cell.b/ny1, cell.c/nz1)
            h_min = min(cell.a/nx1, cell.b/ny1, cell.c/nz1)
            anisotropy = h_max/h_min
            if(anisotropy > 1.1): break
        if i2 == 2:
            st.markdown("reduce grid spacing, anisotropy too large %f"%anisotropy) 
        else:
            i2 = i2//2
            nx1 = (nx+i2-1)//i2 * i2
            ny1 = (ny+i2-1)//i2 * i2
            nz1 = (nz+i2-1)//i2 * i2
        grids_str = col2.text_input("number of grid Nx, Ny, Nz", value="%d %d %d"%(nx1, ny1, nz1))

        hx = cell.a/int(grids_str.split()[0])
        hy = cell.b/int(grids_str.split()[1])
        hz = cell.c/int(grids_str.split()[2])
        st.markdown("final grid spacing: hx =%f hy=%f hz=%f"%(hx,hy,hz))
        anisotropy = max(hx,hy,hz)/min(hx,hy,hz)
        st.markdown("grid anisotropy =%f"%anisotropy)
        pot_grid= col3.text_input("rho pot grid refinement", value="2")

        grid_lines = 'wavefunction_grid="'+grids_str+'"  \n'
        grid_lines += 'potential_grid_refinement="'+pot_grid+'"  \n'

    return grid_lines
