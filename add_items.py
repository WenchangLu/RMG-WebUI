
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
        cstart, col1, col2 = st.columns([0.1,1,1])
        localize_projectors = col1.checkbox("non-local projector localization", value = True, 
                help = "false: non-local projectors will spread in whole space, similar to plane wave codes")
        localize_localpp = col2.checkbox("local potential localization", value = True, 
                help = "fasle: pseudopotential's local part spreads in the whole space")
        max_nlradius = col1.number_input("max radius of non-local projector", 100.0)
        min_nlradius = col2.number_input("max radius of non-local projector", 2.0)
        max_qradius = col1.number_input("max radius of q functions in Ultrasoft PP", 100.0)
        min_qradius = col2.number_input("min radius of q functions in Ultrasoft PP", 2.0)
        write_pseudopotential_plots = col1.checkbox("flag to write pseudopotential plots", False)

        pseudolines += 'localize_localpp ="%s"  \n'%str(localize_localpp)
        pseudolines += 'localize_projectors ="%s"  \n'%str(localize_projectors)
        pseudolines += 'max_nlradius ="%f"  \n'%max_nlradius
        pseudolines += 'min_nlradius ="%f"  \n'%min_nlradius
        pseudolines += 'max_qradius ="%f"  \n'%max_qradius
        pseudolines += 'min_qradius ="%f"  \n'%min_qradius
        pseudolines += 'write_pseudopotential_plots ="%s"  \n'%str(write_pseudopotential_plots)

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
    kp_list_str=col1.text_area("K point list in unit of reciprocal lattice vectors and its weight", "0.0  0.0  0.0  1.0")
    kp_list = kp_list_str.split("\n")
    kpoints = ""
    kpointlines = 'kpoint_mesh = "-1 1 1"  \n'
    num_kpt = 0
    for kp in kp_list:
        if(len(kp.split()) ==4):
            num_kpt+=1
            kpoints += kp + '  \n'
    kpointlines += 'kpoints = "  \n'
    kpointlines += kpoints
    kpointlines += '"  \n'
    col1.markdown(kpointlines)
    if num_kpt == 0:
        st.markdown("kpoint list need to be kx, ky, kz, weight, 4 numbers in a row")
    return kpointlines

def add_kbandstr_lines():
    cs, col1 = st.columns([0.2,1])
    kp_list_str=col1.text_area("special lines for band structure calculation",
       '''   0.0   0.0   0.0   0  G
   0.5   0.0   0.0   20 X''', help ="kx, ky, kz, num, symbol, in unit of reciprocal lattice vector, num: number of kpoinks to previous special kpoint. symbol for plot")
    kpointlines = 'kpoints_bandstructure = "  \n'
    kp_list = kp_list_str.split("\n")
    for kp in kp_list:
        if(len(kp.split()) == 5):
            kpointlines += kp + "  \n"
        else:
            st.markdown("format is wrong for kpoint lines for band structure")

    kpointlines += '"  \n'
    col1.markdown(kpointlines)
    return kpointlines


def add_kpoints(cell):
    expand_ = st.expander("K POINTS")
    with expand_:
        kp_method = st.radio("use gamma point, a mesh or a list", ["gamma", "use mesh", "use list"])
        kp_bandstr = st.radio("kpoints for band structure", ["None", "use special lines", "use list"])
        if kp_method == "gamma":
            kpointlines = 'kpoint_mesh = "1 1 1"  \n'
            kpointlines += 'kpoint_shift = "0 0 0"  \n'
        elif kp_method == "use mesh":
            kpointlines = add_kpoint_mesh(cell)
        if kp_method == "use list" or kp_bandstr == "use list":
            kpointlines = add_kpoint_text()
        if kp_bandstr == "use special lines":    
            kpointlines += add_kbandstr_lines()
    return kpointlines
def add_control():
    expand_ = st.expander("CONTROL OPTIONS")
    extra_lines =""
    with expand_:
        start_mode = st.radio("start mode", 
                ["LCAO Start", "Restart From File", "Random Start",
                  "FIREBALL Start", "Gaussian Start",
                   "Start TDDFT", "Restart TDDFT",
                    "Modified LCAO Start"])
        calculation_mode= st.radio("calculation mode", 
                ["Quench Electrons  ", 
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
        if kohn_sham_solver == "davidson":
            cs, col1,col2,col3 = st.columns([0.1,1,1,1])
            davidson_multiplier = col1.number_input("davidson_multiplier",0)
            davidson_max_steps  = col2.number_input("davidson_max_steps", 8)
            davidson_premg      = col3.number_input("davidson_premg", 4, help = "number of multigrid steps before davidson") 
            extra_lines += 'davidson_multiplier = "%d"  \n'%davidson_multiplier
            extra_lines += 'davidson_max_steps  = "%d"  \n'%davidson_max_steps
            extra_lines += 'davidson_premg      = "%d"  \n'%davidson_premg
        else:
            cs, col1,col2,col3= st.columns([0.1,1,1,1])
            kohn_sham_mg_levels = col1.number_input("kohn_sham_mg_levels", -1, 
                help = "negative: code determines by automatically")
            kohn_sham_pre_smoothing = col2.number_input("kohn_sham_pre_smoothing", 2)
            kohn_sham_post_smoothing = col3.number_input("kohn_sham_post_smoothing", 2)
            kohn_sham_mucycles = col1.number_input("kohn_sham_mucycles", 2)
            kohn_sham_coarse_time_step = col2.number_input("kohn_sham_coarse_time_step", 1.0)
            kohn_sham_time_step = col3.number_input("kohn_sham_time_step", 0.66)
            kohn_sham_mg_timestep = col1.number_input("kohn_sham_mg_timestep", 0.66)

            extra_lines += 'kohn_sham_mg_levels = "%d"  \n'%kohn_sham_mg_levels
            extra_lines += 'kohn_sham_pre_smoothing = "%d"  \n'%kohn_sham_pre_smoothing
            extra_lines += 'kohn_sham_post_smoothing = "%d"  \n'%kohn_sham_post_smoothing
            extra_lines += 'kohn_sham_mucycles = "%d"  \n'%kohn_sham_mucycles
            extra_lines += 'kohn_sham_coarse_time_step = "%f"  \n'%kohn_sham_coarse_time_step
            extra_lines += 'kohn_sham_time_step = "%f"  \n'%kohn_sham_time_step
            extra_lines += 'kohn_sham_mg_timestep = "%f"  \n'%kohn_sham_mg_timestep
        poisson_solver = st.radio("Poisson Solver",
                ["pfft", "multigrid"])
        if poisson_solver == "multigrid":
            cs, col1,col2,col3 = st.columns([0.1,1,1,1])
            poisson_mg_levels = col1.number_input("poisson_mg_levels", -1)
            poisson_pre_smoothing = col2.number_input(" poisson_pre_smoothing", 2)
            poisson_post_smoothing = col3.number_input(" poisson_post_smoothing", 1)
            poisson_mucycles = col1.number_input(" poisson_mucycles", 3)
            poisson_finest_time_step = col2.number_input(" poisson_finest_time_step", 1.0)
            poisson_coarse_time_step = col3.number_input(" poisson_coarse_time_step", 0.8)
            poisson_coarsest_steps = col1.number_input(" poisson_coarsest_steps", 25)
            hartree_max_sweeps = col2.number_input(" hartree_max_sweeps", 10)
            hartree_min_sweeps = col3.number_input(" hartree_min_sweeps", 5)
            extra_lines += 'poisson_mg_levels = "%d"  \n'%poisson_mg_levels
            extra_lines += 'poisson_pre_smoothing = "%d"  \n'% poisson_pre_smoothing
            extra_lines += 'poisson_post_smoothing = "%d"  \n'% poisson_post_smoothing
            extra_lines += 'poisson_mucycles = "%d"  \n'% poisson_mucycles
            extra_lines += 'poisson_finest_time_step = "%f"  \n'% poisson_finest_time_step
            extra_lines += 'poisson_coarse_time_step = "%f"  \n'% poisson_coarse_time_step
            extra_lines += 'poisson_coarsest_steps = "%d"  \n'% poisson_coarsest_steps
            extra_lines += 'hartree_max_sweeps = "%d"  \n'% hartree_max_sweeps
            extra_lines += 'hartree_min_sweeps = "%d"  \n'% hartree_min_sweeps

        subdiag_driver = st.radio("diagonalizatoin libs",
                ["auto", "lapack", "scalapack", "magma", 
                 "cusolver", "elpa", "rocsolver"])

        relax_mass = st.radio("mass for atoms", ["Atomic", "Equal"], 
                help="equal mas for fast relax may help in some cases")
        dos_method = st.radio("density of state calc", 
              ["tetrahedra", "Gaussian"])

        occupations_type = st.radio("occupation type",
                ["Fermi Dirac", "Fixed", "Cold Smearing", "MethfesselPaxton"])
        if occupations_type != "Fixed":
            cs, col1,col2 = st.columns([0.1,1,1])
            occ_smear = col1.number_input("occupation smear in eV", value =0.1)
            MP_order = col2.number_input("Order of Methefessel Paxton Occupation", value=2)


        md_tem_ctrl = st.radio("MD temperature control", 
                ["Nose Hoover Chains","Anderson Rescaling"])
        md_integration_order = st.radio("MD Integration order",
                ["2nd Velocity Verlet",
                 "3rd Beeman-Velocity Verlet",
                 "5th Beeman-Velocity Verlet"])
        md_number_of_nose_thermostats = st.number_input("Number of Nosethermostats", 5)
        ctrl_lines = ""
        ctrl_lines += 'start_mode          ="' +start_mode +'"  \n'
        ctrl_lines += 'calculation_mode    ="' +calculation_mode +'"  \n'
        ctrl_lines += 'kohn_sham_solver    ="' +kohn_sham_solver +'"  \n'
        ctrl_lines += 'subdiag_driver      ="' +subdiag_driver +'"  \n'
        ctrl_lines += 'relax_mass          ="' +relax_mass +'"  \n'
        ctrl_lines += 'dos_method          ="' +dos_method +'"  \n'
        ctrl_lines += 'occupations_type    ="' +occupations_type +'"  \n'
        ctrl_lines += 'occupation_electron_temperature_eV="%f"  \n'%occ_smear
        ctrl_lines += 'MP_order="%d"  \n'%MP_order
        ctrl_lines += 'poisson_solver      ="' +poisson_solver +'"  \n'
        ctrl_lines += 'md_temperature_control    ="' +md_tem_ctrl +'"  \n'
        ctrl_lines += 'md_integration_order="' +md_integration_order +'"  \n'
        ctrl_lines += 'md_number_of_nose_thermostats ="%d"  \n'%md_number_of_nose_thermostats
        ctrl_lines += extra_lines
        return ctrl_lines

def add_grid(cell):
    expand_ = st.expander("REAL SPACE GRID")
    with expand_:
        cs, col1, col2, col3 = st.columns([0.1,1,2,1])
        grid_spacing = col1.number_input("grid spacing(bohr)", value=0.35,
                    help ="use grid spacing to determine the real space grid")
        if cell.unit == "angstrom" :
            grid_spacing = grid_spacing * 0.529177
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
        st.markdown("final grid spacing: hx =%f hy=%f hz=%f "%(hx,hy,hz) + cell.unit)
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
            drho_precond = col2.checkbox("scale q^2/(q^2+q0^2)", False)
            drho_precond_q0= col1.number_input("q0 value", 0.5)
            mixing_lines += 'charge_pulay_Gspace = "' + str(pulay_gspace)+ '"  \n'
            mixing_lines += 'drho_precond = "' +str(drho_precond) + '"  \n'
            mixing_lines += 'drho_precond_q0="%f"  \n'%drho_precond_q0

    return mixing_lines

def add_xc(species):
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
        vdwdf_kernel_filepath = col1.text_input("van der Waals Kernel file", "vdW_kernel_table")
        xc_lines  = 'exchange_correlaton_type="'+xc_type +'"  \n'
        xc_lines += 'exx_mode = "' + exx_mode + '"  \n'
        xc_lines += 'exxdiv_treatment = "' + exxdiv_treatment +'"  \n'
        xc_lines += 'x_gamma_extrapolation ="' + str(x_gamma_extrapolation) +'"  \n'
        xc_lines += 'exx_fracton = "' + exx_fracton +'"  \n'
        xc_lines += 'vdw_corr            ="' +vdw_corr +'"  \n'
        xc_lines += 'vdwdf_grid_type     ="' +vdwdf_grid +'"  \n'
        xc_lines += 'vdwdf_kernel_filepath ="%s"  \n'%vdwdf_kernel_filepath 


        ldaU_mode = st.radio("LDA+U type",["None","Simple"])
        xc_lines += 'ldaU_mode = "%s"  \n'%ldaU_mode
        if(ldaU_mode == "Simple"):
            cs, col1, col2 = st.columns([0.1,1,1])
            Hubbard_U = col1.text_area("HUbbard U for species", "", 
                help= "Ni 6.5 3d 0.0 0.0 0.0 for each specie ")
            xc_lines += 'Hubbard_U ="  \n' + Hubbard_U + '  \n"  \n'
            ldau_mixing_type = col1.radio("mixing type for ldau occupations", ["Linear", "Pulay"])
            xc_lines += 'ldau_mixing_type = "%s"  \n'%ldau_mixing_type
            cs, col1, col2 = st.columns([0.1,1,1])

            ldau_mixing = col1.number_input("mixing fractions", 1.0)
            xc_lines += 'ldau_mixing = "%f"  \n'%ldau_mixing
            if ldau_mixing_type == "Pulay":
                ldau_pulay_order = col2.number_input("Pulay order for lda+u mixing", 5)
                ldau_pulay_scale = col1.number_input("Pulay scale for lda+u mixing", 0.8)
                ldau_pulay_refresh = col2.number_input("Pulay refresh steps for lda+u mixing", 100)
                xc_lines += 'ldau_pulay_order = "%d"  \n'%ldau_pulay_order
                xc_lines += 'ldau_pulay_scale = "%f"  \n'% ldau_pulay_scale
                xc_lines += 'ldau_pulay_refresh = "%d"  \n'%ldau_pulay_refresh

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
            a = sqrt(lattvec[0][0] *lattvec[0][0] +lattvec[0][1] *lattvec[0][1] +lattvec[0][2] *lattvec[0][2] )
            b = sqrt(lattvec[1][0] *lattvec[1][0] +lattvec[1][1] *lattvec[1][1] +lattvec[1][2] *lattvec[1][2] )
            c = sqrt(lattvec[2][0] *lattvec[2][0] +lattvec[2][1] *lattvec[2][1] +lattvec[2][2] *lattvec[2][2] )
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
def add_IOctrl():
    expand_ = st.expander("IO: files and paths")
    with expand_:
        cs, col1, col2 = st.columns([0.2,1,1])
        verbose = col1.checkbox("print out more in log file if True", False)
        cs, col1, col2 = st.columns([0.2,1,1])
        input_wave_function_file = col1.text_input("input wave function file", "Waves/wave.out")
        output_wave_function_file = col2.text_input("output wave function file", "Waves/wave.out")
        write_serial_restart = col1.checkbox("write a serial file for restart", False)
        read_serial_restart = col2.checkbox("restart from a serial file", False)
        compressed_infile = col1.checkbox("read the compressed file for restart", True)
        compressed_outfile = col2.checkbox("compress the out wave file", True)
        input_tddft_file = col1.text_input("input tddft file", "Waves/wave_tddft.out")
        output_tddft_file = col2.text_input("output TDDFT file", "Waves/wave_tddft.out")
        nvme_weights = col1.checkbox("map nonlocal projectors to disk", False)
        nvme_work = col2.checkbox("map work arrays to disk", False)
        nvme_orbitals = col1.checkbox("map orbitals to disk", False)
        nvme_qfunctons = col2.checkbox("map qfunctions to disk", False)
        nvme_weights_filepath = col1.text_input("nvme directory for non-local projectors", "Weights/")
        nvme_work_filepath = col2.text_input("nvme directory for work arrays", "Work/")
        nvme_orbitals_filepath = col1.text_input(" nvme directory for orbitals", "Orbitals/")
        qfunction_filepath = col2.text_input("nvme directory for Qfunction", "Qfunctions/")
        cube_rho = col1.checkbox("output rho in cube format", True)
        output_rho_xsf = col2.checkbox("output rho in xsf format", False)
        cube_vh = col1.checkbox("output vh in cube format",  False)
        cube_pot = col2.checkbox("output pot in cube format", False)
        write_data_period = col1.number_input("steps to write the restart file", 5)
        write_eigvals_period = col2.number_input("steps to write eigenvalues",5)    

    IO_lines = ""
    IO_lines += 'verbose = "%s"  \n'%str(verbose)
    IO_lines += 'input_wave_function_file = "%s"  \n'%input_wave_function_file  
    IO_lines += 'output_wave_function_file = "%s"  \n'%output_wave_function_file 
    IO_lines += 'write_serial_restart = "%s"  \n'%str(write_serial_restart) 
    IO_lines += 'read_serial_restart = "%s"  \n'%str(read_serial_restart) 
    IO_lines += 'compressed_infile = "%s"  \n'%str(compressed_infile) 
    IO_lines += 'compressed_outfile = "%s"  \n'%str(compressed_outfile) 
    IO_lines += 'input_tddft_file = "%s"  \n'%input_tddft_file
    IO_lines += 'output_tddft_file = "%s"  \n'%output_tddft_file
    IO_lines += 'nvme_weights = "%s"  \n'%str(nvme_weights) 
    IO_lines += 'nvme_work = "%s"  \n'%str(nvme_work) 
    IO_lines += 'nvme_orbitals = "%s"  \n'%str(nvme_orbitals) 
    IO_lines += 'nvme_qfunctons = "%s"  \n'%str(nvme_qfunctons) 
    IO_lines += 'nvme_weights_filepath = "%s"  \n'%nvme_weights_filepath
    IO_lines += 'nvme_work_filepath = "%s"  \n'%nvme_work_filepath
    IO_lines += 'nvme_orbitals_filepath = "%s"  \n'%nvme_orbitals_filepath
    IO_lines += 'qfunction_filepath = "%s"  \n'%qfunction_filepath
    IO_lines += 'cube_rho = "%s"  \n'%str(cube_rho) 
    IO_lines += 'output_rho_xsf = "%s"  \n'%str(output_rho_xsf) 
    IO_lines += 'cube_vh = "%s"  \n'%str(cube_vh) 
    IO_lines += 'cube_pot = "%s"  \n'%str(cube_pot) 
    IO_lines += 'write_data_period = "%d"  \n'%write_data_period 
    IO_lines += 'write_eigvals_period = "%d"  \n'%write_eigvals_period
    return IO_lines
def add_spin(species, atoms):
    expand_ = st.expander("SPIN and MAGNETIZATION")
    dict_mag_species = {}
    dict_mag_species = {}
    angle1_species = {}
    angle2_species = {}

    mag = []
    for atom in atoms:
        mag.append([0.0, 0.0, 0.0])
    spin_lines = ""
    with expand_:
        nspin_str = st.radio("spin setup", ["None", "spin polarization", "spin orbit coupling"])

        if(nspin_str == "spin polarization"):
            spin_lines += 'spin_polarization = "True"  \n'
            s_or_a = st.radio("Init Magnetization", ["by species", "by atoms"])
            if s_or_a == "by species":
                for sp in species:
                    dict_mag_species[sp] = st.number_input(sp, 0.0, 
                            help="spin up and down density: (0.5+x, 0.5-x) of total atomic charge density")
                for i in range(len(atoms)):
                    mag[i][0] = dict_mag_species[atoms[i][0]]
            else:
                for i in range(len(atoms)):
                    tem_str = "atom " + str(i) +": up down spin difference"
                    mag[i][0] = st.number_input(tem_str, 0.0)
        
        if(nspin_str == "spin orbit coupling"):
            spin_lines += 'spinorbit = "True"  \n'
            spin_lines += 'noncollinear = "True"  \n'
            s_or_a = st.radio("Init Magnetization", ["by species", "by atoms"])
            cs, col1, col2, col3 = st.columns([0.2,1,1,1])
            if s_or_a == "by species":
                for sp in species:
                    dict_mag_species[sp] = col1.number_input(sp + " mag", 0.0, 
                            help="spin up and down density: (0.5+x, 0.5-x) of total atomic charge density")
                    angle1_species[sp] = col2.number_input(sp + " angle1", 0, help = "180 indicate the -z direction") 
                    angle2_species[sp] = col3.number_input(sp + " angle2", 0, help = "directiopn in xy plane") 
                for i in range(len(atoms)):
                    mag[i][0] = dict_mag_species[atoms[i][0]]
                    mag[i][1] = angle1_species[atoms[i][0]]
                    mag[i][2] = angle2_species[atoms[i][0]]
            else:
                for i in range(len(atoms)):
                    tem_str = "atom " + str(i) +": "
                    mag[i][0] = col1.number_input(tem_str + " mag", 0.0)
                    mag[i][1] = col2.number_input(tem_str + "angle1", 0)
                    mag[i][2] = col3.number_input(tem_str + "angle2", 0)


    return spin_lines, mag



