# Author: Igor Kowalec
# TODO: implement ts_search and vibrate as methods
# TODO: Enable serialization with ASE db
from carmm.run.aims_path import set_aims_command
import os
from ase.io import read


class React_Aims:
    '''Class used to streamline the process of geometry optimisation, input and output generation for
        ASE/FHI-aims setup.'''

    def __init__(self, params: dict, basis_set: str, hpc: str, dimensions: int, filename: str = None):
        '''Define basic parameters'''
        self.data = {}
        self.params = params
        self.hpc = hpc
        self.dimensions = dimensions
        self.basis_set = basis_set
        self.filename = filename

        '''Define additional parameters'''
        self.initial = None
        self.model_optimised = None
        self.model_post_processed = None
        self.final = None
        self.ts = None


    def aims_optimise(self,
                      initial,
                      fmax: float = 0.01,
                      post_process: str = None,
                      relax_unit_cell: bool = False,
                      internal: bool = False,
                      restart: bool = True):

        '''
        TODO: clean up the description
        The function needs information about structure geometry (model), name of hpc system
        to configure FHI-aims environment variables (hpc). Separate directory is created with a naming convention
        based on chemical formula and number of restarts, n (dir_formula_n), ensuring that no outputs are overwritten
        in ASE/FHI-aims.
        The geometry optimisation is restarted from a new Hessian each 80 steps in BFGS algorithm to overcome deep
        potential energy local minima with fmax above convergence criteria. One can choose the type of phase of
        the calculation (gas, surface, bulk) and request a post_processing calculation with a larger basis set.

        PARAMETERS:
        params: dict
            Dictionary containing user's calculator FHI-aims parameters
        initial: Atoms object
        fmax: float
            force convergence criterion for geometry optimisation, i.e. max forces on any atom in eV/A
        post_process: str or None
            Basis set to be used for post_processing if energy calculation using a larger basis set is required
        relax_unit_cell: bool
            True requests a strain filter unit cell relaxation
        restart: bool
            Request restart from previous geometry if True (True by default)

        Returns a list containing the model with data calculated using
        light and tight settings: [model_light, model_tight]
        '''
        from ase.io import read
        from ase.io.trajectory import Trajectory
        from carmm.analyse.forces import is_converged
        from ase.optimize import BFGS

        '''Setup initial parameters'''
        params = self.params
        hpc = self.hpc
        dimensions = self.dimensions
        basis_set = self.basis_set
        self.initial = initial.copy()


        # parent directory
        parent_dir = os.getcwd()

        # Read the geometry
        if self.filename:
            filename = self.filename
        else:
            self.filename = self.initial.get_chemical_formula()

        counter, filename, subdirectory_name = self._restart_setup(self.initial, self.filename, internal, restart)
        out = str(counter) + "_" + str(filename) + ".out"

        # TODO: Perform calculations ONLY if structure is not converged
        # In some instances this could trigger a calculation due to failed ._check_state on the calculator - undesired
        # if not is_converged(model, 0.01):
        if not os.path.exists(subdirectory_name):
            os.mkdir(subdirectory_name)
        os.chdir(subdirectory_name)


        # set the environment variables for geometry optimisation
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020)

        # Restarts will prevent the calculation from getting stuck in deep local minimum when using metaGGA XC mBEEF
        opt_restarts = 0

        if not internal:
            # perform DFT calculations for each filename
            with _calc_generator(params, out_fn=out, dimensions=dimensions, relax_unit_cell=relax_unit_cell)[0] as calculator:
                self.initial.calc = calculator
                while not is_converged(self.initial, fmax):
                    if relax_unit_cell:
                        from ase.constraints import StrainFilter
                        unit_cell_relaxer = StrainFilter(self.initial)
                        opt = BFGS(unit_cell_relaxer,
                                   trajectory=str(counter) + "_" + filename + "_" + str(opt_restarts) + ".traj",
                                   alpha=70.0
                                   )
                    else:
                        opt = BFGS(self.initial,
                                   trajectory=str(counter) + "_" + filename + "_" + str(opt_restarts) + ".traj",
                                   alpha=70.0
                                   )

                    opt.run(fmax=fmax, steps=80)
                    opt_restarts += 1

            os.chdir(parent_dir)

        # setup for internal aims optimiser
        else:
            # perform DFT calculations for each filename
            calculator = _calc_generator(params,
                                 out_fn=out,
                                 dimensions=dimensions,
                                 relax_unit_cell=relax_unit_cell,
                                 internal=internal,
                                 fmax=fmax)

            self.initial.calc = calculator

            if relax_unit_cell:
                # TODO: unit cell relaxtation keyword that works with internal, i.e. relax_unit_cell in FHI-aims
                print("Can't use relax_unit_cell and internal optimisation atm")

            # Just to  trigger energy + force calls
            opt = BFGS(self.initial)
            opt.run(fmax=fmax, steps=0)

            # Save a trajectory from aims output
            from ase.io import Trajectory
            traj = Trajectory(str(counter) + "_" + self.filename + "_" + str(opt_restarts) + ".traj", 'w')
            models = read(out + "@:")

            # TODO: Alert the user of the problem
            # Problem with this approach is inconsistency in energies from aims.out vs communicated over sockets to ASE
            for i in models:
                traj.write(i)

            traj.close()

            # TODO: check if final geometry returned is not identical to initial
            # I think only E+F returned, not changes to structure
            # make sure model returned contains updated geometry
            # model = read(str(counter) + "_" + filename + "_" + str(opt_restarts) + ".traj")

            os.chdir(parent_dir)


        self.model_optimised  = self.initial.copy()

        if post_process:
            print("Commencing calculation using", post_process, "basis set.")

            # Set environment variables for a larger basis set - converged electronic structure
            subdirectory_name_tight = subdirectory_name + "_tight"
            if not os.path.exists(subdirectory_name_tight):
                os.mkdir(subdirectory_name_tight)
            os.chdir(subdirectory_name_tight)
            set_aims_command(hpc=hpc, basis_set=post_process, defaults=2020)

            # Recalculate the structure using a larger basis set in a separate folder
            with _calc_generator(params, out_fn=str(self.filename) + "_" + post_process + ".out",
                         forces=False, dimensions=dimensions)[0] as calculator:
                self.initial.calc = calculator
                self.initial.get_potential_energy()
                traj = Trajectory(self.filename + "_" + post_process + ".traj", "w")
                traj.write(self.initial)
                traj.close()

            # go back to the parent directory to finish the loop
            os.chdir(parent_dir)

            # update the instance with a post_processed model
            self.model_post_processed = self.initial.copy()

            # save the initial geometry also
            self.initial = initial


        return self.model_optimised, self.model_post_processed


    def search_ts(self, initial, final, fmax, unc, interpolation="idpp"):
        '''
        TODO:
        Args:
            initial:
            final:
            fmax:
            unc:
            interpolation:

        Returns:

        '''
        from catlearn.optimize.mlneb import MLNEB

        '''Retrieve common properties'''
        basis_set = self.basis_set
        dimensions = self.dimensions
        hpc = self.hpc
        params = self.params

        # Set the environment parameters
        set_aims_command(hpc=hpc, basis_set=basis_set)

        if not interpolation:
            interpolation = "linear"

        # Read the geometry
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()

        counter, filename, subdirectory_name = self._restart_setup(self.initial, filename)
        out = str(counter) + "_" + str(filename) + ".out"

        if not os.path.exists(subdirectory_name):
            os.mkdir(subdirectory_name)
        os.chdir(subdirectory_name)

        # Desired number of images including start and end point
        # Dense sampling aids convergence but does not increase complexity as significantly as for classic NEB
        n = 0.1

        # Create the sockets calculator - using a with statement means the object is closed at the end.
        with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:
            # Setup the Catlearn object for MLNEB
            neb_catlearn = MLNEB(start=initial,
                                 end=final,
                                 ase_calc=calculator,
                                 n_images=n,
                                 interpolation=interpolation,
                                 mic=True,
                                 restart=True)

            # Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.01 ev/A
            neb_catlearn.run(fmax=fmax,
                             unc_convergence=unc,
                             trajectory='ML-NEB.traj',
                             ml_steps=100,
                             sequential=False,
                             steps=100)


        neb = read("ML-NEB.traj@:")

        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]

        return self.ts


    def _restart_setup(self, model, filename, internal=False, restart=False):
        # TODO: add restart for search_ts
        import fnmatch, shutil
        from ase.io import read

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
        if restart and counter > 0:
            print("Previous optimisation detected", "in dir_" + filename + "_" + str(counter - 1))
            os.chdir(subdirectory_name_prev)

            # check for number of restarted optimisations
            while str(counter - 1) + "_" + filename + "_" + str(opt_restarts) + ".traj" \
                    in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, "*.traj")]:
                opt_restarts += 1

            # read the last optimisation
            if not internal:
                if os.path.exists(str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj"):
                    self.initial = read(str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj")
                    print("Restarting from", str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj")
            else:
                if os.path.exists(str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj"):
                    self.initial = read(str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj")
                    print("Restarting from", str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj")
                elif os.path.exists("geometry.in.next_step"):
                    shutil.copy("geometry.in.next_step", "next_step.in")
                    self.initial = read("next_step.in")
                    print("Restarting from geometry.in.next_step")

            os.chdir("..")

        return counter, filename, subdirectory_name




    def serialize(self):
        '''Save the instance to a file'''
        pass

# TODO: turn into method and include in the class
def _calc_generator(params,
            out_fn="aims.out",
            forces=True,
            dimensions=2,
            sockets=True,
            relax_unit_cell=False,
            internal=False,
            fmax=None):
    '''
    Args:
        out_fn: str
            Alternative name of the FHI-aims output file
        forces: bool
            Flag indicating if force calculation needed
        dimensions: int
            Determines whether we have a "gas"-phase (0) or "periodic" structure (2 or 3)
        sockets: bool
            Flag whether sockets calculator is requested (True by default, way more efficient!). Also lack of sockets
            calculator will result in overwritten .out file as rest of the functionality was written with sockets.
        relax_unit_cell: bool
            Request strain calculation for bulk geometries
        internal: bool
            Using internal FHI-aims optimizer (BFGS-based)
        fmax: float
            force convergence criterion, only in use with internal FHI-aims optimizer
        # TODO: relax_unit_cell and internal should not be used together, Strain Filter is better than unit cell relaxation in FHI-aims
        # TODO: gas does not work with internal due to write_aims producing frac_atoms
    Returns:
        sockets_calc, fhi_calc: sockets calculator and FHI-aims calculator for geometry optimisations
        or
        fhi_calc - FHI-aims calculator for geometry optimisations
    '''

    # New method that gives a default calculator
    import  fnmatch
    if sockets and not internal:
        from carmm.run.aims_calculator import get_aims_and_sockets_calculator
        # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
        # we need to specifically state what the name of the login node is so the two packages can communicate
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=dimensions,
                                                                 verbose=True,
                                                                 codata_warning=False)
    else:
        from carmm.run.aims_calculator import get_aims_calculator
        # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
        # we need to specifically state what the name of the login node is so the two packages can communicate
        fhi_calc = get_aims_calculator(dimensions=dimensions)

    # remove previous xc argument to ensure libxc warning override is first
    fhi_calc.parameters.pop("xc")
    fhi_calc.set(override_warning_libxc='True')

    # Forces required for optimisation:
    if forces:
        fhi_calc.set(compute_forces="true", )

    # add analytical stress keyword for unit cell relaxation
    if relax_unit_cell:
        if not dimensions == 3:
            import sys
            print("Strain Filter calculation requested, but the system is not periodic in all directions.")
            print("If this was intentional, edit geometry_optimisation.py file.")
            print("The calculation will be now terminated.")
            sys.exit()

        fhi_calc.set(compute_analytical_stress='True')

    # use BFGS optimiser internally
    if internal:
        fhi_calc.set(relax_geometry="bfgs " + str(fmax))

    # set a unique .out output name
    fhi_calc.outfilename = out_fn

    # TODO: fundamental settings to be provided as input by the user to make sure nothing essential gets hardcoded
    # FHI-aims settings set up
    fhi_calc.set(**params)

    # For a larger basis set only the Total Potential Energy is required which takes less time than Forces
    if not forces and not internal:
        fhi_calc.set(compute_forces="false",
                     final_forces_cleaned="false")

    if sockets and not internal:
        return sockets_calc, fhi_calc
    else:
        return fhi_calc





"""
def vibrate(params, model, indices, hpc="hawk", dimensions=2, out=None):
    '''
    This method uses ase.vibrations module, see more for info.
    User provides the FHI-aims parameters, the Atoms object and list
    of indices of atoms to be vibrated. FHI-aims environment variables are set
    according to hpc used. Calculation folders are generated automatically
    and a sockets calculator is used for efficiency.

    Work in progress

    Args:
        params: dict
            Dictionary containing user's FHI-aims settings
        model: Atoms object
        indices: list
            List of indices of atoms that require vibrations
        hpc: str
             Name of the hpc facility, valid choices: "hawk", "archer2", "isambard", "young"
        dimensions: int
            0 for gas, 2 for surface and 3 for bulk
        out: str
            Folder names where vibration calculations are kept


    Returns:
        Zero-Point Energy: float

    '''

    from ase.vibrations import Vibrations
    import os, fnmatch

    # set the environment variables for geometry optimisation
    set_aims_command(hpc=hpc, basis_set=basis_set)

    if not out:
        # develop a naming scheme based on chemical formula
        out = model.get_chemical_formula()

    # check/make folders
    counter = 0

    # TODO: check if no aims outputs are overwritten

    # while "vib_" + out + "_" + str(counter) \
    #    in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, "vib_*")]:
    #        counter += 1
    subdirectory_name = "vib_" + out + "_" + str(counter)
    if not os.path.exists(subdirectory_name):
        os.mkdir(subdirectory_name)
    os.chdir(subdirectory_name)

    # calculate vibrations
    with _calc_generator(params, out_fn=out + ".out", dimensions=dimensions)[0] as calculator:
        model.calc = calculator
        vib = Vibrations(model, indices=indices, name=out + "_" + str(counter))
        vib.run()

    # TODO: save trajectories of vibrations visualised for inspection
    vib.summary()
    # vib.write_mode(-1)  # write last mode to trajectory file
    os.chdir("..")

    return vib.get_zero_point_energy()
"""


