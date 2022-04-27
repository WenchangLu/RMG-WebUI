
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
        k_delta = st.number_input("kdelta(2PI/bohr)", value=0.2, help ="use kdelta to estimate kmesh")
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
    with col2:
        kdist = st.text_input("kpoints distribution", value="1", 
                help =" control the parallel over kpoints")
    kpointlines +='kpoint_distribution = "' + kdist +'"  \n'    
    return kpointlines

          

def add_kpoint_text():
    cs, col1 = st.columns([0.2,1])
    kp_list_str=col1.text_area("K point list in unit of reciprocal lattice vectors and its weight")
    kp_list = kp_list_str.split("\n")
    kpoints = ""
    for kp in kp_list:
        if(len(kp.split()) ==4):
            kpoints += kp + '  \n'
    kpointlines = 'kpoints = "  \n'
    kpointlines += kpoints
    kpointlines += '"  \n'
    col1.markdown(kpointlines)
    return kpointlines


def add_kpoints(cell):
    expand_ = st.expander("K POINTS")
    with expand_:
        kp_method = st.radio("use gamma point, a mesh or text input", ["gamma", "use mesh", "text input"])
        if kp_method == "gamma":
            kpointlines = 'kpoint_mesh = "1 1 1"  \n'
            kpointlines += 'kpoint_shift = "0 0 0"  \n'
        elif kp_method == "use mesh":
            kpointlines = add_kpoint_mesh(cell)
        else:
            kpointlines = add_kpoint_text()
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
        ctrl_lines += 'relax_mass          ="' +relax_mass +'"  \n'
        ctrl_lines += 'dos_method          ="' +dos_method +'"  \n'
        ctrl_lines += 'occupations_type    ="' +occupations_type +'"  \n'
        ctrl_lines += 'occupation_electron_temperature_eV="' +occ_smear +'"  \n'
        ctrl_lines += 'poisson_solver      ="' +poisson_solver +'"  \n'
        ctrl_lines += 'md_tem_ctrl         ="' +md_tem_ctrl +'"  \n'
        ctrl_lines += 'md_integration_order="' +md_integration_order +'"  \n'
        return ctrl_lines

def add_grid(cell):
    expand_ = st.expander("REAL SPACE GRID")
    with expand_:
        cs, col1, col2, col3 = st.columns([0.1,1,2,1])
        grid_spacing = col1.number_input("grid spacing(bohr)", value=0.35,
                    help ="use grid spacing to determine the real space grid")
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
        if(anisotropy >=1.1):
            st.markdown('<p style="color:red;">WARNGING: too big grid anisotropy, need to be <1.1 rmg wont run</p>', unsafe_allow_html=True)
        pot_grid= col3.number_input("rho pot grid refinement", value=2)

        grid_lines = 'wavefunction_grid="'+grids_str+'"  \n'
        grid_lines += 'potential_grid_refinement="%d"  \n'%pot_grid

    return grid_lines
def add_scf():
    expand_ = st.expander("SCF & CONVERGENCE CONTROL")
    with expand_:
        cs, col1, col2, col3 = st.columns([0.1,1,1,1])
        max_scf_steps= col1.text_input("max scf steps", value="40")
        max_md_steps= col2.text_input("max md or relax steps", value="10")
        max_exx_steps= col3.text_input("max Exx steps for hybrid or HF", value="20")
        e_err= col1.text_input("energy convergence criterion", value="1.0e-9")
        rms_err= col2.text_input("rms convergence criterion", value="1.0e-7")
        precon_thres= col3.text_input("preconditioner threshold", value="0.0001")
        exx_convergence_criterion = col1.text_input("Exx convergence criterion for hybrid functional", value="1.0e-9") 
        vexx_fft_threshold = col2.text_input("Exx threshold for switch singlt to doulbe precision", value="1.0e-14")
        scf_lines = 'max_scf_steps = "'+max_scf_steps + '"  \n'
        scf_lines += 'max_md_steps = "'+max_md_steps + '"  \n'
        scf_lines += 'max_exx_steps = "'+max_md_steps + '"  \n'
        scf_lines += 'energy_convergence_criterion="' + e_err + '"  \n'
        scf_lines += 'rms_convergence_criterion = "' + rms_err +'"  \n'
        scf_lines += 'preconditioner_threshold = "' + precon_thres + '"  \n'
        scf_lines += 'exx_convergence_criterion = "' + exx_convergence_criterion + '"  \n'
        scf_lines += 'vexx_fft_threshold = "' + vexx_fft_threshold + '"  \n'


    return scf_lines


def add_mixing():
    expand_ = st.expander("MIXING OPTIONS")
    with expand_:
        charge_mixing_type = st.radio("charge density mixing type", 
                ["Broyden", "Pulay", "linear"])
        cs, col1, col2, col3 = st.columns([0.1,1,1,1])
        mix = col1.text_input("charge density mixing parameter",value="0.5") 
        mix_scale = col2.text_input("charge density mixing scale",value="0.5") 
        mix_order = col1.text_input("Broyden or Pulay order", value="5")
        refresh_step = col2.text_input("Broyden or Pulay refresh step", value="100")
        mixing_lines  = 'charge_mixing_type = "'+ charge_mixing_type +'"  \n'
        mixing_lines += 'charge_density_mixing ="' + mix +'"  \n'
        if charge_mixing_type == "Broyden":
            mixing_lines += 'charge_broyden_order = "' + mix_order + '"  \n'
            mixing_lines += 'charge_broyden_scale = "' +mix_scale + '"  \n'
            mixing_lines += 'charge_broyden_refresh = "' +refresh_step + '"  \n'
        elif charge_mixing_type == "Pulay":
            mixing_lines += 'charge_pulay_order = "' + mix_order + '"  \n'
            mixing_lines += 'charge_pulay_scale = "' +mix_scale + '"  \n'
            mixing_lines += 'charge_pulay_refresh = "' +refresh_step + '"  \n'
            pulay_gspace = col1.checkbox("Pulay mixing in G space", False)
            mixing_lines += 'charge_pulay_Gspace = "' + str(pulay_gspace)+ '"  \n'
    return mixing_lines

def add_xc():
    expand_ = st.expander("EXCHANGE CORRELATION POTENTIAL")
    with expand_:
        xc_type = st.radio("exchange correlation type", 
                ["AUTO_XC", "LDA", "GGA XB CP", "PW91", "GGA BLYP", "GGA PBE",
                  "REVPBE", "PW86PBE", "PBESOL", "PBE0", "HSE", "B3LYP", "gaupbe", 
                  "vdw-df", "VDW-DF", "hartree-fock"], 
                help = "AUTO_XC: XC will be determined from pseudopotential")
        cs, col1, col2 = st.columns([0.1,1,1])
        exx_mode = col1.radio("Exx mode", ["Local fft", "Distributed fft"])

        exxdiv_treatment = col2.radio("Exx divergence treatment", 
                ["gygi-baldereschi", "none"])
        x_gamma_extrapolation = col1.checkbox("x_gamma_extrapolation", True)
        exx_fracton = col2.text_input("the fraction of Exx for hybrid functional", value="-1.0", 
                help="negative value: the fraction determined by code for different hybrid functionals")
        vdw_corr = col1.radio("empirical van der Waals correction", 
                ["None", "DFT-D2", "Grimme-D2","DFT-D3"])
        vdwdf_grid = col2.radio("grid for vdw corr",
                ["Coarse", "Fine"])
        xc_lines  = 'exchange_correlaton_type="'+xc_type +'"  \n'
        xc_lines += 'exx_mode = "' + exx_mode + '"  \n'
        xc_lines += 'exxdiv_treatment = "' + exxdiv_treatment +'"  \n'
        xc_lines += 'x_gamma_extrapolation ="' + str(x_gamma_extrapolation) +'"  \n'
        xc_lines += 'exx_fracton = "' + exx_fracton +'"  \n'
        xc_lines += 'vdw_corr            ="' +vdw_corr +'"  \n'
        xc_lines += 'vdwdf_grid_type     ="' +vdwdf_grid +'"  \n'
    return xc_lines

def add_qmcpack():
    expand_ = st.expander("QMCPACK INTERFACE")
    with expand_:
        cs, col1, col2 = st.columns([0.1,1,1])
        qmcpack = col1.checkbox("Write out file for QMCPACK")
        qmcpack_lines = 'write_qmcpack_restart = "' + str(qmcpack) + '"  \n'
        cs, col1, col2 = st.columns([0.1,1,1])
        if qmcpack:
            exx_integrals_filepath = col1.text_input("file name for afqmc", value="afqmc_rmg")
            ExxIntCholosky = col1.checkbox("Cholesky factorization for Vexx", True)
            ExxCholMax = col2.text_input("maximum Cholesky vectors", value="8")
            exx_int_flag = col2.checkbox("Calculate Exack exchange integrals", True)
            qmc_nband = col1.text_input("number of bands for qmcpack", value="0", 
                    help="default value 0: use the number of states")

            qmcpack_lines +='exx_integrals_filepath = "' + exx_integrals_filepath +'"  \n'
            qmcpack_lines +='ExxIntCholosky = "' + str(ExxIntCholosky) +'"  \n'
            qmcpack_lines +='ExxCholMax = "' +  ExxCholMax + '"  \n'
            qmcpack_lines +='exx_int_flag = "' +  str(exx_int_flag) +'"  \n'
            qmcpack_lines +='qmc_nband = "' + qmc_nband +'"  \n'
    return qmcpack_lines


def add_lattice(bounding_box):
    expand_ = st.expander("LATTICE INFO in unit of Anstrom")
    lattvec = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    #estimate the a, b, c = bounding box + 5 Angstrom
    a = bounding_box[1] - bounding_box[0] +5.0
    b = bounding_box[3] - bounding_box[2] +5.0
    c = bounding_box[5] - bounding_box[4] +5.0
    lattvec = [[a,0.0,0.0],[0.0,b,0.0],[0.0,0.0,c]]

    ibrav = 0
    with expand_:
        st.markdown("min_x = %f, max_x = %f"%(bounding_box[0], bounding_box[1]))
        st.markdown("min_y = %f, max_x = %f"%(bounding_box[2], bounding_box[3]))
        st.markdown("min_z = %f, max_z = %f"%(bounding_box[4], bounding_box[5]))
        cs, col1 = st.columns([0.1,1])
        ibrav_str = st.radio("Bravais lattice type", 
                ["Orthorhombic", "Simple Cubic", "FCC", "BCC", "Hexagonal", "do not know"],
                help = "choose do not know for others")
        cs, col1,col2, col3 = st.columns([0.1,1,1,1])
        if ibrav_str == "do not know":
            ibrav = 0
            lattvec_str = col1.text_area("lattice vector in Angstrom", 
                    help = " must be 3x3 numbers")
            mat = lattvec_str.split("\n")
            if len(mat) == 3:
                for i in range(3):
                    vec = mat[i].split()
                    for j in range(3):
                        lattvec[i][j] = float(vec[j])
        elif ibrav_str == "Simple Cubic":
            ibrav = 1
            a = col1.number_input("length a", value=a)
            b = a
            c = a
        elif ibrav_str == "FCC":
            ibrav = 2
            a = col1.number_input("length a", value=a)
            b = a
            c = a
        elif ibrav_str =="BCC":       
            ibrav = 3
            a = col1.number_input("length a", value=a)
            b = a
            c = a
        elif ibrav_str == "Orthorhombic": 
            ibrav = 8
            a = col1.number_input("length a", value=a)
            b = col2.number_input("length b", value=b)
            c = col3.number_input("length c", value=c)
        elif ibrav_str == "Hexagonal":
            ibrav = 4
            a = col1.number_input("length a", value=a)
            b = a
            c = col3.number_input("length c", value=c)
    return (ibrav, a,b,c, lattvec)            
