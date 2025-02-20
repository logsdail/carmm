# Author: Igor Kowalec
import os
#import numpy as np
from ase.calculators.emt import EMT
from ase import Atoms
from ase.io import read
#from ase.vibrations import Vibrations
from carmm.analyse.forces import is_converged
from carmm.run.aims_path import set_aims_command
from ase.io import Trajectory
from carmm.run.workflows.helper import CalculationHelper
from carmm.utils.logger_set import set_logger
from carmm.utils.python_env_check import ase_env_check
#from ase.optimize import BFGS

# TODO: Enable serialization with ASE db - save locations of converged files as well as all properties


class ReactAims:
    """Class used to streamline the process of geometry optimisation, input and output generation for
        ASE/FHI-aims setup."""

    def __init__(self,
                 params: dict,
                 basis_set: str,
                 hpc: str,
                 filename: str = None,
                 nodes_per_instance: int = None,
                 dry_run: bool = False,
                 verbose: bool = True,
                 warning_lvl: int = 1
                 ):
        """
        Args:
            params: dict
                Dictionary containing keywords and parameters for running the FHI-aims calculator
            basis_set: str
                String describing the default basis set as implemented in FHI-aims, e.g. "light" or "tight" etc.
            hpc: str
                Name of the supercomputing facility for setting up environment variables, e.g. "hawk", "archer2"
                See Also run/aims_path
            filename: str
                Naming convention to follow for subsequent calculations and restarts.
                If None the chemical formula is used.
            nodes_per_instance: int
                For parallelised calculation use the above to control the number of nodes for launching each instance
                of FHI-aims.
            dry_run: bool
                A dry run flag for CI-test
            verbose: bool
                Enables the verbosity of the performed operations

        Returns ReactAims object
        """

        """Define basic parameters"""
        self.data = {}
        self.params = params
        self.hpc = hpc
        self.basis_set = basis_set
        self.filename = filename
        self.nodes_per_instance = nodes_per_instance
        self.verbose = verbose

        self.logger = set_logger("react_logger", warning_lvl)

        """Define additional parameters"""
        self.initial = None  # input for optimisation or input for NEB initial image
        self.model_optimised = None  # optimised geometry with calculator attached
        self.model_post_processed = None  # post processed geometry with new calculator attached
        self.final = None  # input final image for NEB
        self.ts = None  # TS geometry from NEB
        self.prev_calcs = None  # for NEB restart
        self.interpolation = None  # TODO: use this for NEBs
        self.dimensions = None

        """ Set the test flag"""
        self.dry_run = dry_run

        if hpc == "custom":
            self.logger.debug(f"WARNING: You have selected 'custom' as an option for HPC.                      ")
            self.logger.debug(f"         This requires a couple of extra steps from the user side.             ")
            self.logger.debug(f"         1) A new environmental variable - CARMM_AIMS_ROOT_DIRECTORY           ")
            self.logger.debug(f"            - must be set. This helps find the folder containing the           ")
            self.logger.debug(f"             default basis function.                                           ")
            self.logger.debug(f"         2) ASE_AIMS_COMMAND must be set, with the correct number of           ")
            self.logger.debug(f"            mpi tasks if desired. Avoid piping output, as React has its own    ")
            self.logger.debug(f"            output folder names.                                               ")

    def aims_optimise(self, atoms: Atoms, fmax: float = 0.01, post_process: str = None, relax_unit_cell: bool = False,
                      restart: bool = True, optimiser=None, opt_kwargs: dict = {}, mace_preopt=False):

        """
         The function needs information about structure geometry (model), name of hpc system
         to configure FHI-aims environment variables (hpc). Separate directory is created with a naming convention
         based on chemical formula and number of restarts, ensuring that no outputs are overwritten
         in ASE/FHI-aims.
         The geometry optimisation is restarted from a new Hessian each 80 steps in BFGS algorithm to overcome deep
         potential energy local minima with fmax above convergence criteria. A post_processing calculation with a
         different basis set can be performed by specifying the basis set name as in post_process.
         PARAMETERS:
         params: dict
             Dictionary containing user's calculator FHI-aims parameters
         atoms: Atoms object
             Contains the geometry information for optimisation
         fmax: float
             Force convergence criterion for geometry optimisation, i.e. max forces on any atom in eV/A
         post_process: str or None
             Basis set to be used for post_processing if energy calculation using a larger basis set is required
         relax_unit_cell: bool
             True requests a strain filter unit cell relaxation
         restart: bool
             Request restart from previous geometry if True (True by default)
         optimiser: optimiser class
            If None - the ase.optimize.BFGS is used
         opt_kwargs: dict
            Dictionary of keyword arguments specific to the provided optimiser class

         Returns a list containing the model with data calculated using your choice of settings
         [model_optimised, model_postprocessed]
         """

        self._initialize_parameters(atoms)

        helper = CalculationHelper(calc_type="Opt",
                                   parent_dir=os.getcwd(),
                                   filename=self.filename,
                                   restart=restart,
                                   verbose=self.verbose)

        counter, out, subdirectory_name, initial = helper.restart_setup()

        if initial is not None:
            self.initial = initial

        if mace_preopt and not initial:
            """Make sure to skip preoptimisations if AIMs restart information is available."""
            assert self._MaceReactor is not None, "Please set MaceReact_Preoptimiser if mace_preopt is True"

            self.initial = self._mace_preoptimise(self.initial, fmax, relax_unit_cell, optimiser, opt_kwargs)

        self._perform_optimization(subdirectory_name, out, counter, fmax, relax_unit_cell, optimiser, opt_kwargs)
        self._finalize_optimization(subdirectory_name, post_process)

        return self.model_optimised, self.model_post_processed

    def _initialize_parameters(self, atoms):
        """
        Internal function for obtaining periodic boundary conditions from the provided Atoms object and generating
         a filename if necessary. Values are then assigned to self.

        Args:
            atoms: Atoms object

        Returns: None
        """
        self.initial = atoms
        self.dimensions = sum(self.initial.pbc)
        if not self.filename:
            self.filename = self.initial.get_chemical_formula()

    def _perform_optimization(self, subdirectory_name: str, out: str, counter: int, fmax: float, relax_unit_cell: bool,
                              optimiser, opt_kwargs: dict):
        """
        An internal function used in aims_optimise to resolve the working directory and perform the optimisation
        calculation of a given structure.

        Args:
            subdirectory_name: string
            out: str
            counter: int
            fmax: float
            relax_unit_cell: bool
            optimiser: bool or optimiser class
            opt_kwargs: dict

        Returns:None

        """
        opt_restarts = 0
        if not optimiser:
            from ase.optimize import BFGS
            optimiser = BFGS

        if not is_converged(self.initial, fmax):
            os.makedirs(subdirectory_name, exist_ok=True)

            """Ensure correct basis set is used"""
            set_aims_command(hpc=self.hpc, basis_set=self.basis_set, defaults=2020,
                             nodes_per_instance=self.nodes_per_instance)

            with _calc_generator(self.params, out_fn=out, dimensions=self.dimensions, relax_unit_cell=relax_unit_cell,
                                    directory=subdirectory_name)[0] as calculator:
                if not self.dry_run:
                    self.initial.calc = calculator
                else:
                    self.initial.calc = EMT()

                while not is_converged(self.initial, fmax):
                    traj_name = f"{subdirectory_name}/{str(counter)}_{self.filename}_{str(opt_restarts)}.traj"

                    if relax_unit_cell:
                        from ase.constraints import StrainFilter
                        unit_cell_relaxer = StrainFilter(self.initial)
                        opt = optimiser(unit_cell_relaxer, trajectory=traj_name, **opt_kwargs)
                    else:
                        opt = optimiser(self.initial, trajectory=traj_name, **opt_kwargs)

                    opt.run(fmax=fmax, steps=80)
                    opt_restarts += 1

            self.model_optimised = read(traj_name)
        else:
            if self.verbose:
                print(f"Structure is converged.")
            self.model_optimised = self.initial

    def _finalize_optimization(self, subdirectory_name: str, post_process: str):
        """
        An internal function used in aims_opimise for postprocessing of an optimisation with a single point energy
        calculation using a different basis set (from the default basis settings in FHI-aims) as specified in the
        post_process variable (e.g. "tight").

        Assigns the post-processed structure to self.model_post_processed.

        Args:
            subdirectory_name: str
            post_process: str

        Returns:
            None
        """

        if post_process:
            subdirectory_name_tight = subdirectory_name + "_" + post_process
            traj_name = f"{subdirectory_name_tight}/{self.filename}_{post_process}.traj"

            if os.path.exists(traj_name):
                self.model_post_processed = read(traj_name)
            else:
                subdirectory_name_tight = subdirectory_name + "_" + post_process
                os.makedirs(subdirectory_name_tight, exist_ok=True)

                set_aims_command(hpc=self.hpc, basis_set=post_process, defaults=2020,
                                 nodes_per_instance=self.nodes_per_instance)

                model_pp = self.model_optimised.copy()

                with _calc_generator(self.params, out_fn=self.filename + "_" + post_process + ".out", forces=False,
                                        dimensions=self.dimensions, directory=subdirectory_name_tight)[0] as calculator:
                    if not self.dry_run:
                        model_pp.calc = calculator
                    else:
                        model_pp.calc = EMT()

                    model_pp.get_potential_energy()
                    traj = Trajectory(traj_name, "w")
                    traj.write(model_pp)
                    traj.close()

                self.model_post_processed = model_pp

    def get_mulliken_charges(self, atoms: Atoms, restart: bool = True):
        """
        This function is used to retrieve atomic charges using Mulliken charge
        decomposition as implemented in FHI-aims. A new trajectory file containing
        the charges is saved.

        Args:
            atoms: Atoms
                Atoms object containing structural information for the calculation
            restart: bool
                If False the calculation is started from scratch, even if previous folders following the naming
                convention are found.

        Returns:
            Atoms object with charges appended
        """

        from ase.io.trajectory import Trajectory
        from carmm.analyse.mulliken import extract_mulliken_charge

        """Setup initial parameters"""
        params = self.params
        hpc = self.hpc
        basis_set = self.basis_set

        self._initialize_parameters(atoms)

        helper = CalculationHelper(calc_type="Charges",
                                   parent_dir=os.getcwd(),
                                   filename=self.filename,
                                   restart=restart,
                                   verbose=self.verbose)

        counter, out, subdirectory_name, initial = helper.restart_setup()

        if initial is not None:
            self.initial = initial
            return self.initial

        """Set the environment variables for geometry optimisation"""
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020)

        """Request Mulliken charge decomposition"""
        params["output"] = ["Mulliken_summary"]

        os.makedirs(subdirectory_name, exist_ok=True)

        with _calc_generator(params, out_fn=out,
                             dimensions=self.dimensions, forces=False,
                             directory=subdirectory_name)[0] as calculator:
            if not self.dry_run:
                self.initial.calc = calculator
            else:
                self.initial.calc = EMT()

            self.initial.get_potential_energy()

        if not self.dry_run:
            charges = extract_mulliken_charge(f'{subdirectory_name}/{out}', len(self.initial))
        else:
            """Return dummy charges for testing"""
            charges = [0 for atom in self.initial]

        self.initial.set_initial_charges(charges)

        traj = Trajectory(f"{subdirectory_name}/{self.filename}_{helper.calc_type.lower()}.traj", 'w')
        traj.write(self.initial)
        traj.close()

        return self.initial

    def search_ts(self, initial: Atoms, final: Atoms,
                  fmax: float, unc: float, interpolation=None,
                  n=0.25, steps=40, restart=True, prev_calcs=None,
                  input_check=0.01, mace_preopt=0, preopt_maxsteps=200):
        """
        This function allows calculation of the transition state using the CatLearn software package in an
        ASE/sockets/FHI-aims setup. The resulting converged band will be located in the MLNEB.traj file.

        Args:
            initial: Atoms object
                Initial structure in the NEB band
            final: Atoms object
                Final structure in the NEB band
            fmax: float
                Convergence criterion of forces in eV/A
            unc: float
                Uncertainty in the fit of the NEB according to the Gaussian Progress Regression model, a secondary
                convergence criterion.
            n: int
                number of middle images, the following is recommended: n * npi = total_no_CPUs
            interpolation: str or []
                The "idpp" or "linear" interpolation types are supported in ASE. alternatively user can provide a custom
                interpolation as a list of Atoms objects.
            n: int or flot
                Desired number of middle images excluding start and end point. If float the number of images is based on
                displacement of atoms. Dense sampling aids convergence but does not increase complexity as significantly
                as for classic NEB.
            restart: bool
                Use previous calculations contained in folders if True, start from scratch if False
            prev_calcs: list of Atoms objects
                Manually provide the training set
            input_check: float or None
                If float the calculators of the input structures will be checked if the structures are below
                the requested fmax and an optimisation will be performed if not.
            mace_preopt: None or str
                Controls whether to use a MACE preoptimised TS path, and the work flow used for preoptimisation
                Requires a MACE calculator be attached to the React_AIMs object via the MaceReact_Preoptimiser
                property.
                None       -     No MACE pre-optimisation
                fullpath   -     MACE preoptimises TS before FHI-aims - FHI-aims inherits all structures from MACE.
                tspath     -     MACE preoptimises after FHI-aims check - MACE receives initial and reactant
                                 structures from FHI-aims and does not optimise
            preopt_maxsteps: int
                Controls the maximum number of steps used in the preoptimiser before giving up and passing the 
                calculation onto MLNEB with FHI-aims.

        Returns: Atoms object
            Transition state geometry structure
        """

        from catlearn.optimize.mlneb import MLNEB
        from ase.io import write

        #TODO: calling mlneb.run() generates files in the current directory, reverting to os.chdir() necessary
        assert not self.nodes_per_instance, "ReactAims.nodes_per_instance is not None \n" \
                                            "Dependency on the catlearn.mlneb module cannot run TS " \
                                            "search concurrently without issues. "

        """Retrieve common properties"""
        basis_set = self.basis_set
        hpc = self.hpc
        params = self.params
        parent_dir = os.getcwd()
        if interpolation:
            self.interpolation = interpolation
        else:
            self.interpolation = "idpp"

        self._initialize_parameters(initial)

        """Set the environment parameters"""
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

        '''Setup the TS calculation'''
        helper = CalculationHelper(calc_type="TS",
                                   parent_dir=os.getcwd(),
                                   filename=self.filename,
                                   restart=restart,
                                   verbose=self.verbose)

        counter, out, subdirectory_name, minimum_energy_path = helper.restart_setup()

        """Set MACE Pre-optimisation flavor"""
        self.mace_preopt_flavour = None
        if mace_preopt is not None and not minimum_energy_path:
            assert isinstance(n, int), "Integer number of images required for MACE TS preoptimiser"
            assert self._MaceReactor is not None, "Please set MaceReact_Preoptimiser if mace_preopt is True."
            assert input_check, "Mace Preoptimisation workflow requires input be set to 'float'."
            self.mace_preopt_flavour = mace_preopt

        """Run preoptimisation flavour 1"""
        if self.mace_preopt_flavour == "fullpath":

            preopt_ts = self._mace_preoptimise_ts(initial, final, fmax, n, self.interpolation,
                                                  input_check, max_steps=preopt_maxsteps)

            initial = preopt_ts[0]
            final = preopt_ts[-1]
            self.mace_interpolation = preopt_ts

        if input_check:
            filename_copy = self.filename
            """Ensure input is converged"""
            if not is_converged(initial, input_check):
                self.filename = f"{filename_copy}_initial"
                initial = self.aims_optimise(initial, input_check, restart=True)[0]
            if not is_converged(final, input_check):
                self.filename =  f"{filename_copy}_final"
                final = self.aims_optimise(final, input_check, restart=True)[0]

            """Set original name after input check is complete"""
            self.filename = filename_copy

        """Run preotimisation flavour 2"""
        if self.mace_preopt_flavour == "tspath":

            preopt_ts = self._mace_preoptimise_ts(initial, final, fmax, n, self.interpolation,
                                                  input_check, max_steps=preopt_maxsteps)

            self.mace_interpolation = preopt_ts

        if not minimum_energy_path:
            minimum_energy_path = [None, None]

        traj_name = f"{subdirectory_name}/ML-NEB.traj"

        if not minimum_energy_path[0]:

            """Create the sockets calculator - using a with statement means the object is closed at the end."""
            with _calc_generator(params,
                                 out_fn=out,
                                 dimensions=self.dimensions,
                                 directory=".")[0] as calculator:  # mlneb files are created in main directory, hence the workaround with current directory and os.chdir

                """Let the user restart from alternative file or Atoms object"""
                if prev_calcs:
                    self.prev_calcs = prev_calcs
                elif minimum_energy_path[1]:
                    self.prev_calcs = minimum_energy_path[1]

                os.makedirs(subdirectory_name, exist_ok=True)

                if self.dry_run:
                    calculator = EMT()

                iterations = 0

                if (self.mace_preopt_flavour is not None) and (not os.path.exists(traj_name)):
                    """GAB: ML-NEB misbehaves if a calculator is not provided for interpolated images"""
                    """     following function ensures correct calculators are attached with closed  """
                    """     sockets.                                                                 """
                    self.interpolation = [ image.copy() for image in self.mace_interpolation[1:-1] ]
                    self.interpolation = [initial] + self.interpolation + [final]

                    for idx, image in enumerate(self.interpolation):
                        self.attach_calculator(self.interpolation[idx], params,
                                           out_fn=out, dimensions=self.dimensions, directory=subdirectory_name,
                                           calc_e=True)

                    initial = self.interpolation[0]
                    final   = self.interpolation[-1]

                while not os.path.exists(traj_name):
                    if iterations > 0:
                        self.prev_calcs = read(f"{subdirectory_name}/last_predicted_path.traj@:")

                    os.chdir(subdirectory_name)

                    """Setup the Catlearn object for MLNEB"""
                    neb_catlearn = MLNEB(start=initial,
                                         end=final,
                                         ase_calc=calculator,
                                         n_images=n,
                                         interpolation=self.interpolation,
                                         neb_method="improvedtangent",
                                         prev_calculations=self.prev_calcs,
                                         mic=True,
                                         restart=restart)
                    if not self.dry_run:
                        """Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.05 eV/A"""
                        neb_catlearn.run(fmax=fmax,
                                         unc_convergence=unc,
                                         trajectory=traj_name,
                                         ml_steps=75,
                                         sequential=False,
                                         steps=steps)

                        iterations += 1
                        os.chdir(parent_dir)
                    else:

                        os.chdir(parent_dir)
                        return None

        """Find maximum energy, i.e. transition state to return it"""
        if minimum_energy_path[0]:
            neb = minimum_energy_path[0]
        else:
            neb = read(f"{traj_name}@:")
        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]

        return self.ts

    def vibrate(self, atoms: Atoms, indices: list, read_only=False):
        """

        This method uses ase.vibrations module, see more for info.
        User provides the FHI-aims parameters, the Atoms object and list
        of indices of atoms to be vibrated. Variables related to FHI-aims are governed by the React object.
        Calculation folders are generated automatically and a sockets calculator is used for efficiency.

        Work in progress

        Args:
            atoms: Atoms object
            indices: list
                List of indices of atoms that require vibrations
            read_only: bool
                Flag for postprocessing - if True, the method only extracts information from existing files,
                no calculations are performed

        Returns:
            Vibrations object
        """
        from ase.vibrations import Vibrations
        import numpy as np

        """Retrieve common properties"""
        basis_set = self.basis_set
        hpc = self.hpc
        params = self.params
        self._initialize_parameters(atoms)
        vib_dir = f"{os.getcwd()}/VibData_{self.filename}/Vibs"
        vib = Vibrations(atoms, indices=indices, name=vib_dir)

        """If a calculation was terminated prematurely (e.g. time limit) empty .json files remain and the calculation
        of the corresponding stretch modes would be skipped on restart. The line below prevents this"""
        vib.clean(empty_files=True)

        """Extract vibration data from existing files"""
        if read_only:
            vib.read()
            return vib
        else:
            """Calculate required vibration modes"""
            required_cache = [os.path.join(vib_dir, "cache." + str(x) + y + ".json") for x in indices for y in [
                "x+", "x-", "y+", "y-", "y-", "z+", "z-"]]
            check_required_modes_files = np.array([os.path.exists(file) for file in required_cache])
            if np.all(check_required_modes_files):
                vib.read()
            else:
                """Set the environment variables for geometry optimisation"""
                set_aims_command(hpc=hpc, basis_set=basis_set,
                                 defaults=2020, nodes_per_instance=self.nodes_per_instance)

                helper = CalculationHelper(calc_type="Vib",
                                           parent_dir=os.getcwd(),
                                           filename=self.filename,
                                           restart=False,
                                           verbose=self.verbose)

                counter, out, subdirectory_name, initial = helper.restart_setup()
                os.makedirs(subdirectory_name, exist_ok=True)

                """Calculate vibrations and write the in a separate directory"""
                with _calc_generator(params, out_fn=out, dimensions=self.dimensions, directory=subdirectory_name
                                     )[0] as calculator:
                    if not self.dry_run:
                        atoms.calc = calculator
                    else:
                        atoms.calc = EMT()
                    vib = Vibrations(atoms, indices=indices, name=vib_dir)
                    vib.run()

            vib.summary()
            vib.write_mode()
            return vib

    @property
    def MaceReact_Preoptimiser(self):
        """
        MACE ASE calculator used in the preoptimisation protocol in aims_optimise or ts_search.

        Args:
             MaceReact MaceReact Obj:
               User-defined MaceReact.

        Returns:
             self._MaceReactor MaceReact Obj:
               User-defined MaceReact.
        """

        return self._MaceReactor

    @MaceReact_Preoptimiser.setter
    def MaceReact_Preoptimiser(self, MaceReact):

        self._MaceReactor = MaceReact

    def _mace_preoptimise(self, atoms: Atoms, fmax, relax_unit_cell, optimiser, opt_kwargs = {}):
        """
        Invokes geometry optimisation method using the MACE ASE calculator defined in self.MaceReact_Preoptimiser.
        The resultant geometry is optimised using FHI-aims in self.aims_optimise. Workflow for optimisation 
        with MACE defined in ReactMACE module.

        Args:
            atoms: Atoms object
            fmax: float
            relax_unit_cell: bool
            optimiser: bool or optimiser class
            opt_kwargs: dict
        Return:
           preopt_atoms: Atoms object
              Atoms object with the calculator object removed
        """

        filname = self._MaceReactor.filename
        self._MaceReactor.filename = "MACE_PREOPT_" + self.filename

        preopt_atoms = self._MaceReactor.mace_optimise(atoms, fmax, restart=False, relax_unit_cell=relax_unit_cell,
                                                       optimiser = optimiser, opt_kwargs = opt_kwargs)

        self._MaceReactor.filename = filname

        # Clean calculator to prevent future problems
        preopt_atoms.calc = None

        return preopt_atoms

    def _mace_preoptimise_ts(self, initial, final, fmax, n, interpolation, input_check, max_steps=200):
        """
        Invokes nudged elastic band (NEB) for an input pathway using the MACE ASE calculator 
        defined in self.MaceReact_Preoptimiser. Workflow for optimisation with MACE defined in
        ReactMACE module.
        
        The resultant reaciton pathway is optimised using FHI-aims in self.search_ts.

        Args:
            initial: Atoms object
            final: Atoms object
            fmax: float
            n: int
            interpolation: string, list of n Atoms objects 
            input_check: None or float

        Returns:
           
        """

        filname = self._MaceReactor.filename
        self._MaceReactor.filename = "MACE_PREOPT_" + self.filename

        if self.mace_preopt_flavour == "fullpath":
            input_check = input_check
        elif self.mace_preopt_flavour == "tspath":
            input_check = None

        preopt_ts = self._MaceReactor.search_ts_neb(initial, final, fmax, n, k=0.05, method="improvedtangent",
                                        interpolation=interpolation, input_check=input_check,
                                        max_steps=max_steps, restart=True)

        self._MaceReactor.filename = filname

        return self._MaceReactor.interpolation.copy()

    def attach_calculator(self, atoms, params,
                    out_fn="aims.out",
                    forces=True,
                    dimensions=2,
                    relax_unit_cell=False,
                    directory=".",
                    calc_e=False):
        """
        Potentially a redundant functionality, but condenses attaching and closing a socket
        calculator to a given Atoms object

        Args:
            calc_e: Boolean, optional
                Calculate total energy through Atoms.get_potential_energy().

        Returns:
            atoms: Atoms
                Input atoms object with correctly closed socket.

        """

        with _calc_generator(params, out_fn=out_fn, forces=forces, dimensions=dimensions,
                             relax_unit_cell=relax_unit_cell,directory=directory)[0] as calculator:

            if self.dry_run:
                calculator = EMT()

            atoms.calc = calculator

            if calc_e:
                atoms.get_potential_energy()

        return atoms

def _calc_generator(params,
                    out_fn="aims.out",
                    forces=True,
                    dimensions=2,
                    relax_unit_cell=False,
                    directory="."):
    """
    This is an internal function for generation of an FHI-aims sockets calculator ensuring that keywords
    required for supported calculation types are added.

    Args:
        out_fn: str
            Alternative name of the FHI-aims output file
        forces: bool
            Flag indicating if force calculation needed
        dimensions: int
            Determines whether we have a "gas"-phase (0) or "periodic" structure (2 or 3)
        relax_unit_cell: bool
            Request strain calculation for bulk geometries
        directory: str
            Set a directory to be used for initialisation of the FHI-aims calculation

    Returns:
        sockets_calc, fhi_calc: sockets calculator and FHI-aims calculator for geometry optimisations
    """

    """New method that gives a default calculator"""

    from carmm.run.aims_calculator import get_aims_and_sockets_calculator
    from ase.calculators.calculator import Parameters
    """On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
    we need to specifically state what the name of the login node is so the two packages can communicate"""
    sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=dimensions,
                                                             logfile=f"{directory}/socketio.log",
                                                             verbose=True,
                                                             codata_warning=False,
                                                             directory=directory,
                                                             **params)

    fhi_calc.parameters['override_warning_libxc'] = 'True'

    """Forces required for optimisation"""
    if not forces:
        fhi_calc.parameters = {k: v for k, v in fhi_calc.parameters.items() if k != 'compute_forces'}

    """Add analytical stress keyword for unit cell relaxation"""
    if relax_unit_cell:
        assert dimensions == 3, "Strain Filter calculation requested, but the system is not periodic in 3 dimensions."
        fhi_calc.parameters['compute_analytical_stress'] = 'True'

    """Sort FHI-aims settings to ensure libxc warning override is prior to xc"""
    keys = list(fhi_calc.parameters.keys())
    keys.sort()
    fhi_calc.parameters = {key: fhi_calc.parameters[key] for key in keys}

    fhi_calc.parameters = Parameters(**fhi_calc.parameters)

    """Set a unique .out output name"""
    if not ase_env_check('3.23.0'):
        fhi_calc.outfilename = out_fn
    else:
        fhi_calc.template.outputname = out_fn

    return sockets_calc, fhi_calc


'''
def serialize(self):
    """Save the instance to a file
    Cases were this would be useful: providing an input file for repeating the calculation using identical settings.
    Such file might be provided alongside code and structures on NOMAD to aid transparency and reproducibility.
    """
    pass

def recover(self):
    """Recover a saved instance form a file"""
'''


