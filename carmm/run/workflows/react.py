# Author: Igor Kowalec
import os
import numpy as np
from ase.calculators.emt import EMT
from ase import Atoms
from ase.io import read
from ase.vibrations import Vibrations
from carmm.analyse.forces import is_converged
from carmm.run.aims_path import set_aims_command
from ase.io import Trajectory
from carmm.run.workflows.helper import CalculationHelper
from ase.optimize import BFGS

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
                 verbose: bool = True):
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
                Enables the the verbosity of the performed operations

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

    def aims_optimise(self, atoms: Atoms, fmax: float = 0.01, post_process: str = None, relax_unit_cell: bool = False,
                      restart: bool = True):

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

        self._perform_optimization(subdirectory_name, out, counter, fmax, relax_unit_cell)
        self._finalize_optimization(subdirectory_name, post_process)

        return self.model_optimised, self.model_post_processed

    def _initialize_parameters(self, atoms):
        """
        Internal function for obtainining periodic boundary conditions from the provided Atoms object and generating
         a filename if necessary. Values are then assigned to self.

        Args:
            atoms: Atoms object

        Returns: None
        """
        self.initial = atoms
        self.dimensions = sum(self.initial.pbc)
        if not self.filename:
            self.filename = self.initial.get_chemical_formula()

    def _perform_optimization(self, subdirectory_name: str, out: str, counter: int, fmax: float, relax_unit_cell: bool):
        """
        An internal function used in aims_optimise to resolve the working directoryand perform the optimisation
        calculation of a given structure.

        Args:
            subdirectory_name: string
            out: str
            counter: int
            fmax: float
            relax_unit_cell: bool

        Returns:None

        """
        opt_restarts = 0

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
                        opt = BFGS(unit_cell_relaxer, trajectory=traj_name, alpha=70.0)
                    else:
                        opt = BFGS(self.initial, trajectory=traj_name, alpha=70.0)

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
            charges = extract_mulliken_charge(out, len(self.initial))
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
                  input_check=0.01):
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

        Returns: Atoms object
            Transition state geometry structure
        """

        from catlearn.optimize.mlneb import MLNEB

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

        '''Setup the TS calculation'''
        helper = CalculationHelper(calc_type="TS",
                                   parent_dir=os.getcwd(),
                                   filename=self.filename,
                                   restart=restart,
                                   verbose=self.verbose)

        counter, out, subdirectory_name, minimum_energy_path = helper.restart_setup()
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
                        return None

        """Find maximum energy, i.e. transition state to return it"""
        if minimum_energy_path[0]:
            neb = minimum_energy_path[0]
        else:
            neb = read(traj_name)
        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]

        return self.ts


    """
    def search_ts_aidneb(self, initial, final, fmax, unc, interpolation=None, n=15,
                         restart=True, prev_calcs=None, input_check=0.01, verbose=True):
        '''
        This function allows calculation of the transition state using the GPAtom software package in an
        ASE/sockets/FHI-aims setup. The resulting converged band will be located in the AIDNEB.traj file.

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
                Desired number of middle images excluding start and end pointpo. If float the number of images is based on
                displacement of atoms. Dense sampling aids convergence but does not increase complexity as significantly
                as for classic NEB.
            restart: bool
                Use previous calculations contained in folders if True, start from scratch if False
            prev_calcs: list of Atoms objects
                Manually provide the training set
            input_check: float or None
                If float the calculators of the input structures will be checked if the structures are below
                the requested fmax and an optimisation will be performed if not.
            verbose: bool
                Flag for turning off printouts in the code

        Returns: Atoms object
            Transition state geometry structure
        '''

        from gpatom.aidneb import AIDNEB

        '''Retrieve common properties'''
        basis_set = self.basis_set
        hpc = self.hpc
        dimensions = sum(initial.pbc)
        params = self.params
        parent_dir = os.getcwd()
        self.interpolation = interpolation

        '''Set the environment parameters'''
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

        if not self.interpolation:
            self.interpolation = "idpp"

        '''Read the geometry'''
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()
            self.filename = filename

        '''Check for previous calculations'''
        counter, subdirectory_name = self._restart_setup("TS", filename, restart=restart, verbose=verbose)

        '''Let the user restart from alternative file or Atoms object'''
        if prev_calcs:
            self.prev_calcs = prev_calcs
            if verbose:
                print("User provided a list of structures manually, training set substituted.")

        elif input_check:
            if not is_converged(initial, input_check):
                self.filename += "_initial"
                initial = self.aims_optimise(initial, input_check, restart=False, verbose=verbose)[0]
                self.initial = self.model_optimised
                '''Set original name after input check is complete'''
                self.filename = filename

            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.aims_optimise(final, input_check, restart=False, verbose=verbose)[0]
                self.final = self.model_optimised
                '''Set original name after input check is complete'''
                self.filename = filename

        out = str(counter) + "_" + str(filename) + ".out"

        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        # TODO: calculating initial and final structure if possible within the GPAtom code

        '''Sockets setup'''
        with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:

            if self.dry_run:
                calculator = EMT()

            '''Training set functionality does not work correctly when reading a list.'''
            '''Instead we save all geometries to a file the GPATOM can use'''
            training_set_dump = Trajectory("AIDNEB_observations.traj", 'w')
            for atoms in self.prev_calcs:
                training_set_dump.write(atoms)
            training_set_dump.close()

            '''Setup the input for AIDNEB'''
            aidneb = AIDNEB(start=initial,
                            end=final,
                            interpolation=self.interpolation,
                            # "idpp" can in some cases (e.g. H2) result in geometry coordinates returned as NaN
                            calculator=calculator,
                            n_images=n+2,
                            max_train_data=40,
                            trainingset="AIDNEB_observations.traj",
                            use_previous_observations=True,
                            neb_method='improvedtangent',
                            mic=True)

            '''Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.01 ev/A'''
            if not self.dry_run:
                aidneb.run(fmax=fmax,
                           unc_convergence=unc,
                           ml_steps=40)
            else:
                os.chdir(parent_dir)
                return None

        '''Find maximum energy, i.e. transition state to return it'''

        neb = read("AIDNEB.traj@" + str(-len(read("initial_path.traj@:")) - 1) + ":") # read last predicted trajectory
        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

        return self.ts

    def search_ts_taskfarm(self, initial, final, fmax, n, method="string", interpolation="idpp", input_check=0.01,
                           max_steps=100, verbose=True):

        '''
        Args:
            initial: Atoms object
                Initial structure in the NEB band
            final: Atoms object
                Final structure in the NEB band
            fmax: float
                Convergence criterion of forces in eV/A
            n: int
                number of middle images, the following is recommended: n * npi = total_no_CPUs
            method: str
                NEB method for the CI-NEB as implemented in ASE, 'string' by default
            interpolation: str or []
                The "idpp" or "linear" interpolation types are supported in ASE. alternatively user can provide a custom
                interpolation as a list of Atoms objects.
            input_check: float or None
                If float the calculators of the input structures will be checked if the structures are below
                the requested fmax and an optimisation will be performed if not.
            max_steps: int
                Maximum number of iteration before stopping the optimizer
            verbose: bool
                Flag for turning off printouts in the code

        Returns: Atoms object
            Transition state geometry structure
        '''

        from ase.neb import NEB
        from ase.optimize import FIRE

        '''Retrieve common properties'''
        basis_set = self.basis_set
        hpc = self.hpc
        dimensions = sum(initial.pbc)
        params = self.params
        parent_dir = os.getcwd()

        '''Set the environment parameters'''
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

        '''Read the geometry'''
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()

        counter, subdirectory_name = self._restart_setup("TS", filename, restart=False, verbose=verbose)

        '''Ensure input is converged'''
        if input_check:
            npi = self.nodes_per_instance
            self.nodes_per_instance = None

            if not is_converged(initial, input_check):
                self.filename = filename + "_initial"
                initial = self.aims_optimise(initial, input_check, restart=False, verbose=False)[0]
            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.aims_optimise(final, input_check, restart=False, verbose=False)[0]

            '''Set original name after input check is complete'''
            self.nodes_per_instance = npi
            self.filename = filename

        out = str(counter) + "_" + str(filename) + ".out"

        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        if interpolation in ["idpp", "linear"]:
            images = [initial]
            for i in range(n):
                image = initial.copy()
                if not self.dry_run:
                    image.calc = _calc_generator(params, out_fn=str(i)+"_"+out, dimensions=dimensions)[0]
                    image.calc.launch_client.calc.directory = "./" + str(i) + "_" + out[:-4]
                else:
                    image.calc = EMT()
                    image.calc.directory = "./" + str(i) + "_" + out[:-4]

                images.append(image)

            images.append(final)
        elif isinstance(interpolation, list):
            assert [isinstance(i, Atoms) for i in interpolation], \
                "Interpolation must be a list of Atoms objects, 'idpp' or 'linear'!"
            assert len(interpolation)-2 == n, \
                "Number of middle images is fed interpolation must match specified n to ensure correct parallelisation"

            images = interpolation
            for i in range(1, len(interpolation)-1):
                # use i-1 for name to retain folder naming as per "idpp"
                if not self.dry_run:
                    images[i].calc = _calc_generator(params, out_fn=str(i-1)+"_"+out, dimensions=dimensions)[0]
                else:
                    images[i].calc = EMT()
                images[i].calc.launch_client.calc.directory = "./"+str(i-1)+"_"+out[:-4]
        else:
            raise ValueError("Interpolation must be a list of Atoms objects, 'idpp' or 'linear'!")

        neb = NEB(images, k=0.05, method=method, climb=True, parallel=True, allow_shared_calculator=False)
        if interpolation in ["idpp", "linear"]:
            neb.interpolate(method=interpolation, mic=True, apply_constraint=True)

        qn = FIRE(neb, trajectory='neb.traj')
        qn.run(fmax=fmax, steps=max_steps)

        for image in images[1:-1]:
            if not self.dry_run:
                image.calc.close()

        '''Find maximum energy, i.e. transition state to return it'''
        self.ts = sorted(images, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

        return self.ts
    """

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


def _calc_generator(params,
                    out_fn="aims.out",
                    forces=True,
                    dimensions=2,
                    relax_unit_cell=False,
                    directory="."):
    """
    This is an internal function for generation of an FHi-aims sockets calculator ensuring that keywords
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
    """On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
    we need to specifically state what the name of the login node is so the two packages can communicate"""
    sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=dimensions,
                                                             logfile=f"{directory}/socketio.log",
                                                             verbose=True,
                                                             codata_warning=False,
                                                             directory=directory
                                                             )

    """Remove previous xc argument to ensure libxc warning override is first"""
    fhi_calc.parameters.pop("xc")
    fhi_calc.set(override_warning_libxc='True')

    """Forces required for optimisation"""
    if not forces:
        fhi_calc.parameters.pop("compute_forces")

    """Add analytical stress keyword for unit cell relaxation"""
    if relax_unit_cell:
        assert dimensions == 3, "Strain Filter calculation requested, but the system is not periodic in all directions."

        fhi_calc.set(compute_analytical_stress='True')

    """Set a unique .out output name"""
    fhi_calc.outfilename = out_fn

    """FHI-aims settings set up"""
    fhi_calc.set(**params)

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


