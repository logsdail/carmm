def my_calc(k_grid,
            out_fn="aims.out",
            forces=True,
            dimensions=2,
            sockets=True,
            preopt=False,
            sf=False,
            internal=False):
    '''
    Args:
        k_grid: (int, int, int)
            Tuple containing 3 positive integers, refers to reciprocal space sampling density
        out_fn: str
            Alternative name of the FHI-aims output file
        forces: bool
            Flag indicating if force calculation needed
        dimensions: int
            Determines whether we have a "gas"-phase (0) or "periodic" structure (2 or 3)
        sockets: bool
            Flag whether sockets calculator is requested (True by default, way more efficient!). Also lack of sockets
            calculator will result in overwritten .out file as rest of the functionality was written with sockets.
        preopt: bool
            Alternative settings choice for methods involving pre-optimisation with a different functional to e.g. find
            metastable geometries easier.
        sf: bool
            Request strain calculation for bulk geometries
        internal: bool
            Using internal FHI-aims optimizer (BFGS-based)
        # TODO: sf and internal should not be used together, Strain Filter is better than unit cell relaxation in FHI-aims
        # TODO: gas does not work with internal due to write_aims producing frac_atoms
    Returns:
        sockets_calc, fhi_calc: sockets calculator and FHI-aims calculator for geometry optimisations
        or
        fhi_calc - FHI-aims calculator for geometry optimisations

    '''
    # New method that gives a default calculator
    import os, fnmatch
    if sockets and not internal:
        from carmm.run.aims_calculator import get_aims_and_sockets_calculator
        # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
        # we need to specifically state what the name of the login node is so the two packages can communicate
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=dimensions,
                                                                 k_grid=k_grid,
                                                                 verbose=True,
                                                                 codata_warning=False)
    else:
        from carmm.run.aims_calculator import get_aims_calculator
        # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
        # we need to specifically state what the name of the login node is so the two packages can communicate
        fhi_calc = get_aims_calculator(dimensions=dimensions,
                                       k_grid=k_grid)

    # remove previous xc argument to ensure libxc warning override is first
    fhi_calc.parameters.pop("xc")
    fhi_calc.set(override_warning_libxc='True')

    # FHI-aims settings set up
    fhi_calc.set(xc_pre=['pbe', '10'],
                 xc='libxc MGGA_X_MBEEF+GGA_C_PBE_SOL',
                 spin='none',
                 relativistic=('atomic_zora','scalar'),
                 compute_forces="true",
                 #output=['mulliken'],
                 #use_dipole_correction='True',
                 #force_correction="true",
                 sc_accuracy_etot=1e-6,
                 sc_accuracy_eev=1e-3,
                 sc_accuracy_rho=1e-6,
                 sc_accuracy_forces=1e-4,
                 final_forces_cleaned='true',
                 )

    # Use full PBEsol exchange-correlation to preptimise geometries
    if preopt:
        fhi_calc.set(xc='libxc GGA_X_PBE_SOL+GGA_C_PBE_SOL')

    # add analytical stress keyword for unit cell relaxation
    if sf:
        if not dimensions == 3:
            import sys
            print("Strain Filter calculation requested, but the system is not periodic in all directions.")
            print("If this was intentional, edit geometry_optimisation.py file.")
            print("The calculation will be now terminated.")
            sys.exit()

        fhi_calc.set(compute_analytical_stress = 'True')

    # use BFGS optimiser internally
    if internal:
        fhi_calc.set(relax_geometry="bfgs 0.010")


    # For a larger basis set only the Total Potential Energy is required which takes less time than Forces
    if not forces and not internal:
        fhi_calc.set(compute_forces="false",
                     final_forces_cleaned="false")


    #counter = 0
    # check/make folders
    #while str(counter) + "_" + out_fn \
    #        in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, out_fn)]:
    #    counter += 1
    #fhi_calc.outfilename = str(counter) + "_" + out_fn

    # in case one would want multiple .out files in one directory, though that makes life difficult for output parsers
    # e.g. in NOMAD repostiory

    fhi_calc.outfilename = out_fn

    if sockets and not internal:
        return sockets_calc, fhi_calc
    else:
        return fhi_calc


def aims_optimise(model,
                  hpc,
                  constraints=False,
                  dimensions=2,
                  fmax=0.01,
                  tight=True,
                  preopt=False,
                  sf=False,
                  internal=False):
    '''
    Function used to streamline the process of geometry optimisation, input and ouput generation for ASE/FHI-aims setup.
    The function needs information about structure geometry (model), name of hpc system to configure FHI-aims
    environment variables (hpc). Separate directory is created with a naming convention based on chemical
    formula and number of restarts, n (dir_formula_n), ensuring that no outputs are overwritten in ASE/FHI-aims.
    The geometry optimisation is restarted after each 10 descent steps in BFGSLineSearch algorithm to overcome deep
    potential energy local minima which result in fmax above convergence criteria. One can choose the type of phase of
    the calculation (gas, surface, bulk) and request a postprocessing calculation with a larger basis set.

    PARAMETERS:
    model: Atoms object
    hpc: string
        Name of the hpc facility, valid choices: "hawk", "archer2", "isambard", "young"
    constraints:
        List of constraint objects - see ase.constraints (e.g. FixAtoms)
    dimensions: int
        0, 2, 3 for gas, surface and bulk calculations, respectively
        See Also carmm.run.get_aims_calculator
    fmax: float
        force convergence criterion for geometry optimisation, i.e. max forces on any atom in eV/A
    tight: bool
        Set to True if postprocessing energy calculation using a larger basis set is required
    preopt: bool
        If true, an alternative set of settings will be passed to the FHI-aims calculator
        See Also my_calc()
    sf: bool
        True requests a strain filter unit cell relaxation

    Returns a list containing the model with data calculated using
    light and tight settings: [model_light, model_tight]
    '''

    from carmm.run.aims_path import set_aims_command
    from ase.io import read
    import os, fnmatch, shutil
    from ase.io.trajectory import Trajectory
    from carmm.analyse.forces import is_converged
    from ase.optimize import BFGSLineSearch

    # Read the geometry
    filename = model.get_chemical_formula()

    # define k_grid
    k_grid = get_k_grid(model, 0.018, dimensions=dimensions, verbose=True)

    # ensure separate folders are in place for each calculation input
    counter = 0
    # Count .traj outputs
    opt_restarts = 0

    # check/make folders
    while "dir_" + filename + "_" + str(counter) \
            in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, "dir_*")]:
        
        counter += 1

    subdirectory_name = "dir_" + filename + "_" + str(counter)
    subdirectory_name_prev = "dir_" + filename + "_" + str(counter - 1)

    # check previous calculations for convergence
    if counter > 0:
        print("Previous optimisation detected", "in dir_" + filename + "_" + str(counter - 1))
        os.chdir(subdirectory_name_prev)

        # check for number of restarted optimisations
        while str(counter - 1) + "_" + filename + "_" + str(opt_restarts) + ".traj" \
                in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, "*.traj")]:
            opt_restarts += 1

        # read the last optimisation
        if not internal:
            if os.path.exists(str(counter - 1) + "_" + filename + "_" + str(opt_restarts-1) + ".traj"):
                model = read(str(counter - 1) + "_" + filename + "_" + str(opt_restarts-1) + ".traj")
                print("Restarting from", str(counter - 1) + "_" + filename + "_" + str(opt_restarts-1) + ".traj")
        else:
            if os.path.exists(str(counter-1) + "_" + filename + "_" + str(opt_restarts-1) + ".traj"):
                model = read(str(counter-1) + "_" + filename + "_" + str(opt_restarts-1) + ".traj")
                print("Restarting from", str(counter-1) + "_" + filename + "_" + str(opt_restarts-1) + ".traj")
            elif os.path.exists("geometry.in.next_step"):
                shutil.copy("geometry.in.next_step", "next_step.in")
                model = read("next_step.in")
                print("Restarting from geometry.in.next_step")

        os.chdir("..")

    # TODO: Perform calculations ONLY if structure is not converged
    # if not is_converged(model, 0.01):
    if not os.path.exists(subdirectory_name):
        os.mkdir(subdirectory_name)
    os.chdir(subdirectory_name)

    out = str(counter) + "_" + str(filename) + ".out"

    # set the environment variables for geometry optimisation
    set_aims_command(hpc=hpc, basis_set="light")

    # Restart optimisations after 10 descent steps in BFGSLineSearch
    # Repeat until converged (10 descent steps is usually around 30-50 force evals)
    # Restarts will prevent the calculation from getting stuck in deep local minimum when using metaGGA XC mBEEF
    opt_restarts = 0

    if not internal:
        # perform DFT calculations for each filename
        with my_calc(k_grid, out_fn=out, dimensions=dimensions, preopt=preopt, sf=sf)[0] as calculator:
            model.calc = calculator
            if constraints:
                model.set_constraint(constraints)

            if sf:
                from ase.constraints import StrainFilter
                model = StrainFilter(model)

            while not is_converged(model, fmax):
                opt = BFGSLineSearch(model,
                           trajectory=str(counter) + "_" + filename + "_" + str(opt_restarts) + ".traj",
                           maxstep=0.2,
                           alpha=50.0
                            )

                opt.run(fmax=fmax, steps=10)
                opt_restarts += 1

        os.chdir("..")

    # setup for internal aims optimiser
    else:
        # perform DFT calculations for each filename
        calculator = my_calc(k_grid,
                 out_fn=out,
                 dimensions=dimensions,
                 preopt=preopt,
                 sf=sf,
                 internal=internal)

        model.calc = calculator

        if constraints:
            model.set_constraint(constraints)

        if sf:
            print("Can't use sf and internal optimisation atm")

        # Just to  trigger energy + force calls
        opt = BFGSLineSearch(model)
        opt.run(fmax=fmax, steps=0)

        # Save a trajectory from aims output
        from ase.io import Trajectory
        traj = Trajectory(str(counter) + "_" + filename + "_" + str(opt_restarts) + ".traj", 'w')
        models = read(out+"@:")

        for i in models:
            traj.write(i)

        traj.close()
        # make sure model returned contains updated geometry
        #model = read(str(counter) + "_" + filename + "_" + str(opt_restarts) + ".traj")



        os.chdir("..")

    if tight:
        print("Commencing calculation using 'tight' basis set.")
        # store a copy of the converged model with "light" calculator data
        model_tight = model.copy()

        # Set environment variables for a larger basis set - converged electronic structure
        subdirectory_name_tight = subdirectory_name + "_tight"
        if not os.path.exists(subdirectory_name_tight):
            os.mkdir(subdirectory_name_tight)
        os.chdir(subdirectory_name_tight)
        set_aims_command(hpc=hpc, basis_set="tight")

        # Recalculate the structure using a larger basis set in a separate folder
        with my_calc(k_grid, out_fn=str(filename) + "_tight.out",
                    forces=False, dimensions=dimensions)[0] as calculator:
            model_tight.calc = calculator
            model_tight.get_potential_energy()
            traj = Trajectory(filename + "_tight.traj", "w")
            traj.write(model_tight)
            traj.close()

        # go back to the parent directory to finish the loop
        os.chdir("..")

        return [model, model_tight]
    else:
        return [model]



def get_k_grid(model, sampling_density, dimensions=2, verbose=False):

    '''
    Based converged value of reciprocal space sampling provided,
    this function analyses the xyz-dimensions of the simulation cell
    and returns the minimum k-grid as a tuple that can be used
    as a value for the k_grid keyword in the Aims calculator.
    This function allows to maintain consistency in input
    for variable supercell sizes.

    Parameters:
    model: Atoms object
        Periodic model that requires k-grid for calculation in FHI-aims.
    sampling_density: float
        Converged value of minimum reciprocal space sampling required for
        accuracy of the periodic calculation. Value is a fraction between
        0 and 1, unit is /Å.
    dimensions: int
        2 sets the k-grid in z-direction to 1 for surface slabs, 3 calculates as normal, k_grid not necessary for others
        that have vacuum padding added.
    verbose: bool
        Flag turning print statements on/off

    Returns:
        float containing 3 integers: (kx, ky, kz)
        or
        None if a non-periodic model is presented
    '''


    import math
    import numpy as np
    # define k_grid sampling density /A
    x = np.linalg.norm(model.get_cell()[0])
    y = np.linalg.norm(model.get_cell()[1])
    z = np.linalg.norm(model.get_cell()[2])

    if np.all([x, y, z] == 0):
        return False

    k_x = math.ceil((1 / sampling_density) * (1 / x))
    k_y = math.ceil((1 / sampling_density) * (1 / y))
    # recognise surface models and set k_z to 1
    if dimensions == 2:
        k_z = 1
    elif dimensions == 3:
        k_z = math.ceil((1 / sampling_density) * (1 / z))
    else:
        print("Wrong number of periodic dimensions specified.")
        return False

    k_grid = (k_x, k_y, k_z)

    if verbose:
        print("Based on lattice xyz dimensions", "x", round(x, 3), "y", round(y, 3), "z", round(z, 3))
        print("and", str(sampling_density), "sampling density, the k-grid chosen for periodic calculation is",
              str(k_grid) + ".")
        print()

    return k_grid


def vibrate(model, indices, hpc="hawk", dimensions=2, out=None):
    """Class for calculating vibrational modes using finite difference.

    The vibrational modes are calculated from a finite difference
    approximation of the Hessian matrix.

    The *summary()*, *get_energies()* and *get_frequencies()* methods all take
    an optional *method* keyword.  Use method='Frederiksen' to use the method
    described in:

      T. Frederiksen, M. Paulsson, M. Brandbyge, A. P. Jauho:
      "Inelastic transport theory from first-principles: methodology and
      applications for nanoscale devices", Phys. Rev. B 75, 205413 (2007)

    atoms: Atoms object
        The atoms to work on.
    indices: list of int
        List of indices of atoms to vibrate.  Default behavior is
        to vibrate all atoms.
    name: str
        Name to use for files.
    delta: float
        Magnitude of displacements.
    nfree: int
        Number of displacements per atom and cartesian coordinate, 2 and 4 are
        supported. Default is 2 which will displace each atom +delta and
        -delta for each cartesian coordinate.
    """

    from ase.vibrations import Vibrations
    from carmm.run.aims_path import set_aims_command
    import os, fnmatch

    # set the environment variables for geometry optimisation
    set_aims_command(hpc=hpc, basis_set="light")

    if not out:
        # develop a naming scheme based on chemical formula
        out = model.get_chemical_formula()

    # check/make folders
    counter = 0
    # while "vib_" + out + "_" + str(counter) \
    #    in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, "vib_*")]:
    #        counter += 1
    subdirectory_name = "vib_" + out + "_" + str(counter)
    if not os.path.exists(subdirectory_name):
        os.mkdir(subdirectory_name)
    os.chdir(subdirectory_name)

    # calculate vibrations
    with my_calc(get_k_grid(model, 0.018, dimensions=dimensions, verbose=True), out_fn=out + ".out", dimensions=dimensions)[0] as calculator:
        model.calc = calculator
        vib = Vibrations(model, indices=indices, name=out + "_" + str(counter))
        vib.run()

    # TODO: save trajectories of vibrations visualised
    vib.summary()
    # vib.write_mode(-1)  # write last mode to trajectory file
    os.chdir("..")

    return vib.get_zero_point_energy()
