# Author: Igor Kowalec
import fnmatch
import os
import numpy as np

from ase.calculators.emt import EMT
from ase import Atoms
from ase.io import read
from ase.vibrations import Vibrations
from carmm.analyse.forces import is_converged
from carmm.run.aims_path import set_aims_command


# TODO: Enable serialization with ASE db - save locations of converged files as well as all properties
# TODO: rework the use of filename with the calculator to use calc.directory instead as recommended by ASE authors


class ReactAims:
    """Class used to streamline the process of geometry optimisation, input and output generation for
        ASE/FHI-aims setup."""

    def __init__(self,
                 params: dict,
                 basis_set: str,
                 hpc: str,
                 filename: str = None,
                 nodes_per_instance: int = None,
                 dry_run: bool = False):
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

        Returns ReactAims object
        """

        """Define basic parameters"""
        self.data = {}
        self.params = params
        self.hpc = hpc
        self.basis_set = basis_set
        self.filename = filename
        self.nodes_per_instance = nodes_per_instance

        """Define additional parameters"""
        self.initial = None                 # input for optimisation or input for NEB initial image
        self.model_optimised = None         # optimised geometry with calculator attached
        self.model_post_processed = None    # post processed geometry with new calculator attached
        self.final = None                   # input final image for NEB
        self.ts = None                      # TS geometry from NEB
        self.prev_calcs = None              # for NEB restart

        """ Set the test flag"""
        self.dry_run = dry_run

    def aims_optimise(self,
                      atoms: Atoms,
                      fmax: float = 0.01,
                      post_process: str = None,
                      relax_unit_cell: bool = False,
                      restart: bool = True,
                      verbose: bool = True):

        """
        The function needs information about structure geometry (model), name of hpc system
        to configure FHI-aims environment variables (hpc). Separate directory is created with a naming convention
        based on chemical formula and number of restarts, n (opt_formula_n), ensuring that no outputs are overwritten
        in ASE/FHI-aims.
        The geometry optimisation is restarted from a new Hessian each 80 steps in BFGS algorithm to overcome deep
        potential energy local minima with fmax above convergence criteria. One can choose the type of phase of
        the calculation (gas, surface, bulk) and request a post_processing calculation with a larger basis set.

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

        Returns a list containing the model with data calculated using
        light and tight settings: [model_light, model_tight]
        """
        from ase.io import read
        from ase.io.trajectory import Trajectory
        from carmm.analyse.forces import is_converged
        from ase.optimize import BFGS

        """Setup initial parameters"""
        params = self.params
        hpc = self.hpc
        basis_set = self.basis_set
        self.initial = atoms
        dimensions = sum(self.initial.pbc)
        i_geo = atoms.copy()
        i_geo.calc = atoms.calc

        """Parent directory"""
        parent_dir = os.getcwd()

        """Read the geometry"""
        if not self.filename:
            self.filename = self.initial.get_chemical_formula()

        filename = self.filename

        counter, subdirectory_name = self._restart_setup("Opt", self.filename, restart, verbose=verbose)
        out = str(counter) + "_" + str(filename) + ".out"

        """Perform calculation only if required"""
        if is_converged(atoms, fmax):
            if verbose:
                print("The forces are below", fmax, "eV/A. No calculation required.")
            self.model_optimised = self.initial
        elif is_converged(self.initial, fmax):
            if verbose:
                print("The forces are below", fmax, "eV/A. No calculation required.")
            self.model_optimised = self.initial
            self.initial = i_geo
        else:
            os.makedirs(subdirectory_name, exist_ok=True)
            os.chdir(subdirectory_name)

            """Set the environment variables for geometry optimisation"""
            set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

            """Occasional optimizer restarts will prevent the calculation from getting stuck in deep local minimum"""
            opt_restarts = 0

            """Perform DFT calculations for each filename"""
            with _calc_generator(params,
                                 out_fn=out,
                                 dimensions=dimensions,
                                 relax_unit_cell=relax_unit_cell)[0] as calculator:

                if not self.dry_run:
                    self.initial.calc = calculator
                else:
                    self.initial.calc = EMT()

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

            self.model_optimised = read(str(counter) + "_" + filename + "_" + str(opt_restarts-1) + ".traj")
            os.chdir(parent_dir)

        self.initial = i_geo

        if post_process:
            if verbose:
                print("Commencing calculation using", post_process, "basis set.")

            model_pp = self.model_optimised.copy()

            """Set environment variables for a larger basis set - converged electronic structure"""
            subdirectory_name_tight = subdirectory_name + "_" + post_process
            os.makedirs(subdirectory_name_tight, exist_ok=True)
            os.chdir(subdirectory_name_tight)

            set_aims_command(hpc=hpc, basis_set=post_process, defaults=2020, nodes_per_instance=self.nodes_per_instance)

            """Recalculate the structure using a larger basis set in a separate folder"""
            with _calc_generator(params, out_fn=str(self.filename) + "_" + post_process + ".out",
                                 forces=False, dimensions=dimensions)[0] as calculator:
                if not self.dry_run:
                    model_pp.calc = calculator
                else:
                    model_pp.calc = EMT()
                    
                model_pp.get_potential_energy()
                traj = Trajectory(self.filename + "_" + post_process + ".traj", "w")
                traj.write(model_pp)
                traj.close()

            """Go back to the parent directory to finish the loop"""
            os.chdir(parent_dir)

            """update the instance with a post_processed model"""
            self.model_post_processed = model_pp

        return self.model_optimised, self.model_post_processed

    def get_mulliken_charges(self, initial: Atoms, verbose=True):

        """
        This function is used to retrieve atomic charges using Mulliken charge
        decomposition as implemented in FHI-aims. A new trajectory file containing
        the charges

        Args:
            initial: Atoms
                Atoms object containing structural information for the calculation
            verbose: bool
                Flag for turning off printouts in the code

        Returns:
            Atoms object with charges appended
        """

        from ase.io.trajectory import Trajectory
        from carmm.analyse.mulliken import extract_mulliken_charge

        """Setup initial parameters"""
        params = self.params
        hpc = self.hpc
        basis_set = self.basis_set
        self.initial = initial
        dimensions = sum(self.initial.pbc)

        """Parent directory"""
        parent_dir = os.getcwd()

        """Read the geometry"""
        if not self.filename:
            self.filename = self.initial.get_chemical_formula()

        filename = self.filename
        assert type(filename) == str, "Invalid type, filename should be string"

        counter, subdirectory_name = self._restart_setup("Charges", self.filename)

        """Check for previously completed calculation"""
        if os.path.exists(os.path.join(subdirectory_name[:-1]+str(counter-1), filename+"_charges.traj")):
            file_location = os.path.join(subdirectory_name[:-1]+str(counter-1), filename+"_charges.traj")
            self.initial = read(file_location)
            if verbose:
                print("Previously calculated structure has been found at", file_location)
            return self.initial

        out = str(counter) + "_" + str(filename) + ".out"

        """Set the environment variables for geometry optimisation"""
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020)

        """Request Mulliken charge decomposition"""
        params["output"] = ["Mulliken_summary"]

        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        with _calc_generator(params, out_fn=out, dimensions=dimensions, forces=False)[0] as calculator:
            if not self.dry_run:
                self.initial.calc = calculator
            else:
                self.initial.calc = EMT()

            self.initial.get_potential_energy()

        if not self.dry_run:
            charges = extract_mulliken_charge(out, len(self.initial))
        else:
            charges = initial.get_charges()

        self.initial.set_initial_charges(charges)

        traj = Trajectory(filename+"_charges.traj", 'w')
        traj.write(self.initial)
        traj.close()

        os.chdir(parent_dir)

        return self.initial

    def search_ts(self, initial, final,
                  fmax, unc, interpolation="idpp",
                  n=0.25, restart=True, prev_calcs=None,
                  input_check=0.01, verbose=True):
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
            verbose: bool
                Flag for turning off printouts in the code

        Returns: Atoms object
            Transition state geometry structure
        """

        from catlearn.optimize.mlneb import MLNEB

        """Retrieve common properties"""
        basis_set = self.basis_set
        hpc = self.hpc
        dimensions = sum(initial.pbc)
        params = self.params
        parent_dir = os.getcwd()

        """Set the environment parameters"""
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

        if not interpolation:
            interpolation = "idpp"

        """Read the geometry"""
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()
            self.filename = filename

        counter, subdirectory_name = self._restart_setup("TS", filename, restart=restart, verbose=verbose)
        if os.path.exists(os.path.join(subdirectory_name[:-1] + str(counter-1), "ML-NEB.traj")):
            previously_converged_ts_search = os.path.join(subdirectory_name[:-1] + str(counter-1), "ML-NEB.traj")
            print("TS search already converged at", previously_converged_ts_search)

            neb = read(previously_converged_ts_search+"@:")
            self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]
            os.chdir(parent_dir)

            return self.ts

        elif input_check:
            """Ensure input is converged"""
            if not is_converged(initial, input_check):
                self.filename += "_initial"
                initial = self.aims_optimise(initial, input_check, restart=True, verbose=False)[0]
            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.aims_optimise(final, input_check, restart=True, verbose=False)[0]

            """Set original name after input check is complete"""
            self.filename = filename

        out = str(counter) + "_" + str(filename) + ".out"

        """Let the user restart from alternative file or Atoms object"""
        if prev_calcs:
            self.prev_calcs = prev_calcs

        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        """Create the sockets calculator - using a with statement means the object is closed at the end."""
        with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:
            if self.dry_run:
                calculator = EMT()
            iterations = 0
            while not os.path.exists('ML-NEB.traj'):
                if iterations > 0:
                    self.prev_calcs = read("last_predicted_path.traj@:")
                    interpolation = self.prev_calcs

                """Setup the Catlearn object for MLNEB"""
                neb_catlearn = MLNEB(start=initial,
                                     end=final,
                                     ase_calc=calculator,
                                     n_images=n,
                                     interpolation=interpolation,
                                     neb_method="improvedtangent",
                                     prev_calculations=self.prev_calcs,
                                     mic=True,
                                     restart=restart)
                if not self.dry_run:
                    """Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.05 eV/A"""
                    neb_catlearn.run(fmax=fmax,
                                     unc_convergence=unc,
                                     trajectory='ML-NEB.traj',
                                     ml_steps=75,
                                     sequential=False,
                                     steps=40)

                    iterations += 1
                else:
                    os.chdir(parent_dir)
                    return None

        """Find maximum energy, i.e. transition state to return it"""
        neb = read("ML-NEB.traj@:")
        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

        return self.ts

    def search_ts_aidneb(self, initial, final, fmax, unc, interpolation="idpp", n=15,
                         restart=True, prev_calcs=None, input_check=0.01, verbose=True):
        """
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
            verbose: bool
                Flag for turning off printouts in the code

        Returns: Atoms object
            Transition state geometry structure

        """
        from gpatom.aidneb import AIDNEB

        """Retrieve common properties"""
        basis_set = self.basis_set
        hpc = self.hpc
        dimensions = sum(initial.pbc)
        params = self.params
        parent_dir = os.getcwd()

        """Set the environment parameters"""
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

        if not interpolation:
            interpolation = "idpp"

        """Read the geometry"""
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()
            self.filename = filename

        """Check for previous calculations"""
        counter, subdirectory_name = self._restart_setup("TS", filename, restart=restart, verbose=verbose)

        """Let the user restart from alternative file or Atoms object"""
        if prev_calcs:
            self.prev_calcs = prev_calcs
            if verbose:
                print("User provided a list of structures manually, training set substituted.")

        elif input_check:
            if not is_converged(initial, input_check):
                self.filename += "_initial"
                initial = self.aims_optimise(initial, input_check, restart=False, verbose=verbose)[0]
                self.initial = self.model_optimised
                """Set original name after input check is complete"""
                self.filename = filename

            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.aims_optimise(final, input_check, restart=False, verbose=verbose)[0]
                self.final = self.model_optimised
                """Set original name after input check is complete"""
                self.filename = filename

        out = str(counter) + "_" + str(filename) + ".out"

        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        # TODO: calculating initial and final structure if possible within the GPAtom code

        """Sockets setup"""
        with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:

            if self.dry_run:
                calculator = EMT()

            """Setup the input for AIDNEB"""
            aidneb = AIDNEB(start=initial,
                            end=final,
                            interpolation=interpolation,
                            # "idpp" can in some cases (e.g. H2) result in geometry coordinates returned as NaN
                            calculator=calculator,
                            n_images=n+2,
                            max_train_data=50,
                            trainingset=self.prev_calcs,
                            use_previous_observations=True,
                            neb_method='improvedtangent',
                            mic=True)

            """Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.01 ev/A"""
            if not self.dry_run:
                aidneb.run(fmax=fmax,
                           unc_convergence=unc,
                           ml_steps=100)
            else:
                os.chdir(parent_dir)
                return None

        """Find maximum energy, i.e. transition state to return it"""
        neb = read("AIDNEB.traj@:")
        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

        return self.ts

    def search_ts_taskfarm(self, initial, final, fmax, n, method="string", interpolation="idpp", input_check=0.01,
                           max_steps=100, verbose=True):
        """

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

        """
        from ase.neb import NEB
        from ase.optimize import FIRE

        """Retrieve common properties"""
        basis_set = self.basis_set
        hpc = self.hpc
        dimensions = sum(initial.pbc)
        params = self.params
        parent_dir = os.getcwd()

        """Set the environment parameters"""
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

        """Read the geometry"""
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()

        counter, subdirectory_name = self._restart_setup("TS", filename, restart=False, verbose=verbose)

        """Ensure input is converged"""
        if input_check:
            npi = self.nodes_per_instance
            self.nodes_per_instance = None

            if not is_converged(initial, input_check):
                self.filename = filename + "_initial"
                initial = self.aims_optimise(initial, input_check, restart=False, verbose=False)[0]
            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.aims_optimise(final, input_check, restart=False, verbose=False)[0]

            """Set original name after input check is complete"""
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

        """Find maximum energy, i.e. transition state to return it"""
        self.ts = sorted(images, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

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
            Zero-Point Energy: float
        """

        """Retrieve common properties"""
        basis_set = self.basis_set
        hpc = self.hpc
        params = self.params
        parent_dir = os.getcwd()
        dimensions = sum(atoms.pbc)

        if not self.filename:
            """develop a naming scheme based on chemical formula"""
            self.filename = atoms.get_chemical_formula()

        vib_dir = parent_dir + "/VibData_" + self.filename + "/Vibs"
        print(vib_dir)

        vib = Vibrations(atoms, indices=indices, name=vib_dir)

        """If a calculation was terminated prematurely (e.g. time limit) empty .json files remain and the calculation
        of the corresponding stretch modes would be skipped on restart. The line below prevents this"""
        vib.clean(empty_files=True)

        """Extract vibration data from existing files"""
        if read_only:
            vib.read()

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

                """Generate a unique folder for aims calculation"""
                counter, subdirectory_name = self._restart_setup("Vib",
                                                                 filename=self.filename,
                                                                 restart=False,
                                                                 verbose=False)

                os.makedirs(subdirectory_name, exist_ok=True)
                os.chdir(subdirectory_name)

                """Name the aims output file"""
                out = str(counter) + "_" + str(self.filename) + ".out"

                """Calculate vibrations and write the in a separate directory"""
                with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:
                    if not self.dry_run:
                        atoms.calc = calculator
                    else:
                        atoms.calc = EMT()

                    vib = Vibrations(atoms, indices=indices, name=vib_dir)
                    vib.run()

            vib.summary()

        """Generate a unique folder for aims calculation"""
        if not read_only:
            os.chdir(vib_dir)
            vib.write_mode()
        os.chdir(parent_dir)

        return vib.get_zero_point_energy()

    def _restart_setup(self, calc_type, filename, restart=False, verbose=True):
        """This is an internal function for generation of an FHi-aims working folders and ensuring that DFT outputs are
        not overwritten.

        Args:
            calc_type: str
                "Opt", "Vib", "TS", "Charges" are currently supported choices
                TODO: separate TS search types into separate restarts
            filename: str
                String containing the naming scheme of the files from calculation runs
            restart: bool
                If True, files from previous calculation runs will be searched and read for calculation restart
            verbose: bool
                Flag for turning code printouts on and off

        Returns:
            counter: int
                Number of current run if previous calculations are detected
            subdirectory_name: str
                Name of the current update calculation run


        """
        _supported_calc_types = ["Opt", "Vib", "TS", "Charges"]
        assert calc_type in _supported_calc_types
        calc_type += "_"

        """Ensure separate folders are in place for each calculation input"""
        counter = 0

        """Check/make folders"""
        parent_dir = os.getcwd()

        while calc_type + filename + "_" + str(counter) \
                in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, calc_type + "*")]:
            counter += 1

        folder_counter = counter
        subdirectory_name = calc_type + filename + "_" + str(counter)
        subdirectory_name_prev = calc_type + filename + "_" + str(counter - 1)

        """Check previous calculations for convergence"""
        if restart and counter > 0:
            restart_found = False
            while not restart_found and counter > 0:


                if verbose:
                    print("Previous calculation detected in", calc_type + filename + "_" + str(counter - 1))
                os.chdir(calc_type + filename + "_" + str(counter - 1))

                if calc_type == "TS_":
                    if os.path.exists("evaluated_structures.traj"):
                        self.prev_calcs = read("evaluated_structures.traj@:")
                        restart_found = True
                        break
                    elif os.path.exists("AID_observations.traj"):
                        self.prev_calcs = read("AID_observations.traj@:")
                        restart_found = True
                        break
                    elif verbose:
                        print('Previous trajectory not found, starting from scratch.')

                elif calc_type == "Opt_":
                    """Check for number of restarted optimisations"""
                    opt_restarts = 0
                    while os.path.exists(str(counter - 1) + "_" + filename + "_" + str(opt_restarts) + ".traj"):
                        opt_restarts += 1

                    """Read the last optimisation"""
                    traj_name = str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj"
                    while os.path.exists(traj_name):
                        if os.path.getsize(traj_name):
                            self.initial = read(traj_name)
                            restart_found = True
                            if verbose:
                                print("Restarting from", traj_name)
                            break

                        elif verbose:
                            print(traj_name + ".traj file empty!")

                        traj_name = str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj"
                        opt_restarts -= 1

                counter -= 1
                os.chdir("..")

        os.chdir(parent_dir)

        return folder_counter, subdirectory_name

    '''
    def serialize(self):
        """Save the instance to a file"""
        pass

    def recover(self):
        """Recover a saved instance form a file"""
    '''


# TODO: turn into method and include in the class
def _calc_generator(params,
                    out_fn="aims.out",
                    forces=True,
                    dimensions=2,
                    relax_unit_cell=False):
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

    Returns:
        sockets_calc, fhi_calc: sockets calculator and FHI-aims calculator for geometry optimisations
    """

    """New method that gives a default calculator"""

    from carmm.run.aims_calculator import get_aims_and_sockets_calculator
    """On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
    we need to specifically state what the name of the login node is so the two packages can communicate"""
    sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=dimensions,
                                                             verbose=True,
                                                             codata_warning=False
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
