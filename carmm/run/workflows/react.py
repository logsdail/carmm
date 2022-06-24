# Author: Igor Kowalec
import fnmatch
import os
import shutil

import numpy as np
from ase import Atoms
from ase.io import read
from ase.vibrations import Vibrations
from carmm.analyse.forces import is_converged
from carmm.run.aims_path import set_aims_command


# TODO: Enable serialization with ASE db - save locations of converged files as well as all properties
# TODO: rework the use of filename with the calculator to use calc.directory instead as recommended by ASE authors


class ReactAims:
    '''Class used to streamline the process of geometry optimisation, input and output generation for
        ASE/FHI-aims setup.'''

    def __init__(self,
                 params: dict,
                 basis_set: str,
                 hpc: str,
                 filename: str = None,
                 nodes_per_instance: int = None):

        '''Define basic parameters'''
        self.data = {}
        self.params = params
        self.hpc = hpc
        self.basis_set = basis_set
        self.filename = filename
        self.nodes_per_instance = nodes_per_instance # n nodes used enabling parallel calcs

        '''Define additional parameters'''
        self.initial = None # input for optimisation or input for ML-NEB initial image
        self.model_optimised = None # optimised geometry with calculator attached
        self.model_post_processed = None # post processed geometry with new calculator attached
        self.final = None # input final image for ML-NEB
        self.ts = None # TS geometry from ML-NEB
        self.prev_calcs = None # for ML-NEB restart


    def aims_optimise(self,
                      atoms: Atoms,
                      fmax: float = 0.01,
                      post_process: str = None,
                      relax_unit_cell: bool = False,
                      internal: bool = False,
                      restart: bool = True,
                      verbose: bool = True):

        '''
        TODO: clean up the description
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
        basis_set = self.basis_set
        self.initial = atoms
        dimensions = sum(self.initial.pbc)
        i_geo = atoms.copy()
        i_geo.calc = atoms.calc



        # parent directory
        parent_dir = os.getcwd()

        # Read the geometry
        if not self.filename:
            self.filename = self.initial.get_chemical_formula()

        filename = self.filename

        counter, subdirectory_name = self._restart_setup("Opt", self.filename, internal, restart, verbose=verbose)
        out = str(counter) + "_" + str(filename) + ".out"


        '''Perform calculation only if required'''
        if is_converged(atoms, fmax):
            print("The forces are below", fmax, "eV/A. No calculation required.")
            self.model_optimised = self.atoms
        elif is_converged(self.initial, fmax):
            print("The forces are below", fmax, "eV/A. No calculation required.")
            self.model_optimised = self.initial
            self.initial = i_geo
        else:
            os.makedirs(subdirectory_name, exist_ok=True)
            os.chdir(subdirectory_name)


            # set the environment variables for geometry optimisation
            set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

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

                self.model_optimised = read(str(counter) + "_" + filename + "_" + str(opt_restarts-1) + ".traj")
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

                if relax_unit_cell and verbose:
                    # TODO: unit cell relaxation keyword that works with internal, i.e. relax_unit_cell in FHI-aims
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
                self.model_optimised = read(str(counter) + "_" + self.filename + "_" + str(opt_restarts) + ".traj")

                # TODO: check if final geometry returned is not identical to initial
                # I think only E+F returned, not changes to structure
                # make sure model returned contains updated geometry
                # model = read(str(counter) + "_" + filename + "_" + str(opt_restarts) + ".traj")

                os.chdir(parent_dir)


        self.initial = i_geo

        if post_process:
            if verbose:
                print("Commencing calculation using", post_process, "basis set.")

            model_pp = self.model_optimised.copy()

            # Set environment variables for a larger basis set - converged electronic structure
            subdirectory_name_tight = subdirectory_name + "_" + post_process
            os.makedirs(subdirectory_name_tight, exist_ok=True)
            os.chdir(subdirectory_name_tight)

            set_aims_command(hpc=hpc, basis_set=post_process, defaults=2020, nodes_per_instance=self.nodes_per_instance)

            # Recalculate the structure using a larger basis set in a separate folder
            with _calc_generator(params, out_fn=str(self.filename) + "_" + post_process + ".out",
                         forces=False, dimensions=dimensions)[0] as calculator:
                model_pp.calc = calculator
                model_pp.get_potential_energy()
                traj = Trajectory(self.filename + "_" + post_process + ".traj", "w")
                traj.write(model_pp)
                traj.close()

            # go back to the parent directory to finish the loop
            os.chdir(parent_dir)

            # update the instance with a post_processed model
            self.model_post_processed = model_pp

        return self.model_optimised, self.model_post_processed


    def get_mulliken_charges(self, initial: Atoms, verbose=True):

        '''
        This function is used to retrieve atomic charges using Mulliken charge
        decomposition as implemented in FHI-aims. A new trajectory file containing
        the charges

        Args:
            initial: Atoms
                Atoms object containing structural information for the calculation

        Returns:
            Atoms object with charges appended

        '''

        from ase.io.trajectory import Trajectory
        from carmm.analyse.mulliken import extract_mulliken_charge

        '''Setup initial parameters'''
        params = self.params
        hpc = self.hpc
        basis_set = self.basis_set
        self.initial = initial
        dimensions = sum(self.initial.pbc)

        # parent directory
        parent_dir = os.getcwd()

        # Read the geometry
        if not self.filename:
            self.filename = self.initial.get_chemical_formula()

        filename = self.filename
        assert type(filename) == str, "Invalid type, filename should be string"

        counter, subdirectory_name = self._restart_setup("Charges", self.filename)

        # Check for previously completed calculation
        if os.path.exists(os.path.join(subdirectory_name[:-1]+str(counter-1), filename+"_charges.traj")):
            file_location = os.path.join(subdirectory_name[:-1]+str(counter-1), filename+"_charges.traj")
            self.initial = read(file_location)
            if verbose:
                print("Previously calculated structure has been found at", file_location)
            return self.initial

        out = str(counter) + "_" + str(filename) + ".out"

        # set the environment variables for geometry optimisation
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020)

        # request Mulliken charge decomposition
        params["output"]= ["Mulliken_summary"]


        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        with _calc_generator(params, out_fn=out, dimensions=dimensions, forces=False)[0] as calculator:
            self.initial.calc = calculator
            self.initial.get_potential_energy()

        charges = extract_mulliken_charge(out, len(self.initial))
        # charges = [round(charge, 2) for charge in charges] # rounding not necessary in original file
        self.initial.set_initial_charges(charges)

        traj = Trajectory(filename+"_charges.traj", 'w')
        traj.write(self.initial)
        traj.close()

        os.chdir(parent_dir)

        return self.initial


    def search_ts(self, initial, final, fmax, unc, interpolation="idpp", restart=True, prev_calcs=None, input_check=0.01, verbose=True):
        '''

        Args:
            initial: Atoms object
            final: Atoms object
            fmax: float
            unc: float
            interpolation: str or list of Atoms objects
            restart: bool
            prev_calcs: list of Atoms objects
            input_check: float or None
            verbose: bool

        Returns: Atoms object
            Transition state geometry structure

        '''
        from catlearn.optimize.mlneb import MLNEB

        '''Retrieve common properties'''
        basis_set = self.basis_set
        hpc = self.hpc
        dimensions = sum(initial.pbc)
        params = self.params
        parent_dir = os.getcwd()

        # Set the environment parameters
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

        if not interpolation:
            interpolation = "idpp"

        # Read the geometry
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

        # ensure input is converged
        elif input_check:

            if not is_converged(initial, input_check):
                self.filename += "_initial"
                initial = self.aims_optimise(initial, input_check, restart=False, verbose=False)[0]
            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.aims_optimise(final, input_check, restart=False, verbose=False)[0]

            # Set original name after input check is complete
            self.filename = filename

        out = str(counter) + "_" + str(filename) + ".out"

        # Let the user restart from alternative file or Atoms object
        if prev_calcs:
            self.prev_calcs = prev_calcs

        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        # Desired number of images including start and end point
        # Dense sampling aids convergence but does not increase complexity as significantly as for classic NEB
        n = 0.25

        # Create the sockets calculator - using a with statement means the object is closed at the end.
        with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:
            iterations = 0
            while not os.path.exists('ML-NEB.traj'):
                if iterations > 0:
                    self.prev_calcs = read("last_predicted_path.traj@:")
                    interpolation = self.prev_calcs

                # Setup the Catlearn object for MLNEB
                neb_catlearn = MLNEB(start=initial,
                                     end=final,
                                     ase_calc=calculator,
                                     n_images=n,
                                     #k=0.05,
                                     interpolation=interpolation,
                                     neb_method="aseneb",
                                     prev_calculations=self.prev_calcs,
                                     mic=True,
                                     restart=restart)
                # TODO: TEST

                # Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.01 ev/A
                neb_catlearn.run(fmax=fmax,
                                 unc_convergence=unc,
                                 trajectory='ML-NEB.traj',
                                 ml_steps=75,
                                 sequential=False,
                                 steps=40)

                iterations += 1

        # Find maximum energy, i.e. transition state to return it
        neb = read("ML-NEB.traj@:")
        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

        return self.ts


    def search_ts_aidneb(self, initial, final, fmax, unc, interpolation="idpp",
                         restart=True, prev_calcs=None, input_check=0.01, verbose=True):
        '''

        Args:
            initial: Atoms object
            final: Atoms object
            fmax: float
            unc: float
            interpolation: str or list of Atoms objects
            restart: bool
            prev_calcs: list of Atoms objects
            input_check: float or None
            verbose: bool

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

        # Set the environment parameters
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

        if not interpolation:
            interpolation = "idpp"

        # Read the geometry
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()
            self.filename = filename

        counter, subdirectory_name = self._restart_setup("TS", filename, restart=restart, verbose=verbose)
        if os.path.exists(os.path.join(subdirectory_name[:-1] + str(counter-1), "AIDNEB.traj")):
            previously_converged_ts_search = os.path.join(subdirectory_name[:-1] + str(counter-1), "AIDNEB.traj")
            print("TS search already converged at", previously_converged_ts_search)

            neb = read(previously_converged_ts_search+"@:")
            self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]
            os.chdir(parent_dir)

            return self.ts

        # ensure input is converged
        elif input_check:

            if not is_converged(initial, input_check):
                self.filename += "_initial"
                initial = self.aims_optimise(initial, input_check, restart=False, verbose=False)[0]
                self.initial = self.model_optimised
            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.aims_optimise(final, input_check, restart=False, verbose=False)[0]
                self.final = self.model_optimised

            # Set original name after input check is complete
            self.filename = filename

        out = str(counter) + "_" + str(filename) + ".out"

        # Let the user restart from alternative file or Atoms object
        if prev_calcs:
            self.prev_calcs = prev_calcs

        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        # Desired number of images including start and end point
        # Dense sampling aids convergence but does not increase complexity as significantly as for classic NEB
        n = 0.25

        # TODO: calculating initial and final structure if possible within the GPAatom code

        '''sockets setup'''
        with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:

            # Setup the GPAtom object for AIDNEB
            aidneb = AIDNEB(start=initial,
                            end=final,
                            interpolation=interpolation,
                            # "idpp" can in some cases (e.g. H2) result in geometry coordinates returned as NaN, no error exit, but calculator stuck
                            calculator=calculator,
                            n_images=n,
                            max_train_data=50,
                            trainingset=self.prev_calcs,
                            use_previous_observations=True,
                            neb_method='improvedtangent',
                            mic=True)

            # Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.01 ev/A
            aidneb.run(fmax=fmax,
                       unc_convergence=unc,
                       ml_steps=100)


        '''non-sockets setup
        calculator = _calc_generator(params, out_fn=out, dimensions=dimensions, sockets=False)

        # Setup the GPAtom object for AIDNEB
        aidneb = AIDNEB(start=initial,
                        end=final,
                        interpolation=interpolation, # "idpp" can in some cases (e.g. H2) result in geometry coordinates returned as NaN, no error, but calculator stuck
                        calculator=calculator,
                        n_images=n,
                        trainingset=self.prev_calcs,
                        use_previous_observations=True,
                        mic=True) #
        # TODO: TEST

        # Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.01 ev/A
        aidneb.run(fmax=fmax,
                    unc_convergence=unc,
                    ml_steps=100)
        '''

        # Find maximum energy, i.e. transition state to return it
        neb = read("AIDNEB.traj@:")
        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

        return self.ts



    def search_ts_taskfarm(self, initial, final, fmax, n, method="string", interpolation="idpp", input_check=0.01, verbose=True):
        '''

        Args:
            initial: Atoms object
            final: Atoms object
            fmax: float
            n: int
                number of images, the following is recommended: n * npi = total_no_CPUs
            interpolation: str
            input_check: float or None
            verbose: bool

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

        # Set the environment parameters
        set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)


        # Read the geometry
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()

        counter, subdirectory_name = self._restart_setup("TS", filename, restart=False, verbose=verbose)

        # ensure input is converged
        if input_check:
            npi = self.nodes_per_instance
            self.nodes_per_instance = None

            if not is_converged(initial, input_check):
                self.filename = filename + "_initial"
                initial = self.aims_optimise(initial, input_check, restart=False, verbose=False)[0]
            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.aims_optimise(final, input_check, restart=False, verbose=False)[0]

            # Set original name after input check is complete
            self.nodes_per_instance = npi
            self.filename = filename

        out = str(counter) + "_" + str(filename) + ".out"

        os.makedirs(subdirectory_name, exist_ok=True)
        os.chdir(subdirectory_name)

        if interpolation in ["idpp", "linear"]:
            images = [initial]
            for i in range(n):
                image = initial.copy()
                image.calc = _calc_generator(params, out_fn=str(i)+"_"+out, dimensions=dimensions)[0]
                image.calc.launch_client.calc.directory = "./"+str(i)+"_"+out[:-4]
                images.append(image)

            images.append(final)
        elif isinstance(interpolation, list):
            assert [isinstance(i, Atoms) for i in interpolation], "Interpolation must be a list of Atoms objects, 'idpp' or 'linear'!"
            assert len(interpolation)-2 == n, "Number of middle images is fed interpolation must match specified n to ensure correct parallelisation"

            images = interpolation
            for i in range(1, len(interpolation)-1):
                # use i-1 for name to retain folder naming as per "idpp"
                images[i].calc = _calc_generator(params, out_fn=str(i-1)+"_"+out, dimensions=dimensions)[0]
                images[i].calc.launch_client.calc.directory = "./"+str(i-1)+"_"+out[:-4]
        else:
            raise ValueError("Interpolation must be a list of Atoms objects, 'idpp' or 'linear'!")


        neb = NEB(images, k=0.05, method=method, climb=True, parallel=True, allow_shared_calculator=False)
        if interpolation in ["idpp", "linear"]:
            neb.interpolate(method=interpolation, mic=True, apply_constraint=True)

        qn = FIRE(neb, trajectory='neb.traj')
        qn.run(fmax=fmax, steps=100)

        for image in images[1:-1]:
            image.calc.close()

        # Find maximum energy, i.e. transition state to return it
        self.ts = sorted(images, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

        return self.ts


    def vibrate(self, atoms: Atoms, indices: list, read_only=False):
        '''

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
        '''

        '''Retrieve common properties'''
        basis_set = self.basis_set
        hpc = self.hpc
        params = self.params
        parent_dir = os.getcwd()
        dimensions = sum(atoms.pbc)

        if not self.filename:
            # develop a naming scheme based on chemical formula
            self.filename = atoms.get_chemical_formula()

        vib_dir = parent_dir + "/VibData_" + self.filename +"/Vibs"
        print(vib_dir)

        vib = Vibrations(atoms, indices=indices, name=vib_dir)
        # If a calculation was terminated prematurely (e.g. time limit) empty .json files remain and the calculation
        # of the corresponding stretch modes would be skipped on restart. The line below prevents this
        vib.clean(empty_files=True)
        # Extract vibration data from existing files
        if read_only:
            vib.read()
        # Calculate required vibration modes
        else:
            required_cache = [os.path.join(vib_dir,"cache."+str(x)+y+".json") for x in indices for y in [
                "x+" ,"x-", "y+", "y-", "y-", "z+", "z-"]]
            check_required_modes_files = np.array([os.path.exists(file) for file in required_cache])

            if np.all(check_required_modes_files == True):
                vib.read()
            else:
                # set the environment variables for geometry optimisation
                set_aims_command(hpc=hpc, basis_set=basis_set, defaults=2020, nodes_per_instance=self.nodes_per_instance)

                # check/make folders
                counter = 0
                # Generate a unique folder for aims calculation
                counter, subdirectory_name = self._restart_setup("Vib",
                                                                 filename=self.filename,
                                                                 restart=False,
                                                                 verbose=False)

                os.makedirs(subdirectory_name, exist_ok=True)
                os.chdir(subdirectory_name)

                # Name the aims output file
                out = str(counter) + "_" + str(self.filename) + ".out"

                # calculate vibrations and write the in a separate directory
                with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:
                    atoms.calc = calculator
                    vib = Vibrations(atoms, indices=indices, name=vib_dir)
                    vib.run()

            vib.summary()


        # Generate a unique folder for aims calculation
        if not read_only:
            os.chdir(vib_dir)
            vib.write_mode()
        os.chdir(parent_dir)

        return vib.get_zero_point_energy()


    def _restart_setup(self, calc_type, filename, internal=False, restart=False, verbose=True):
        '''

        Args:
            calc_type:
            filename:
            internal:
            restart:
            verbose:

        Returns:

        '''
        # TODO: add restart for search_ts_dyneb
        _supported_calc_types = ["Opt", "Vib", "TS", "Charges"]
        assert calc_type in _supported_calc_types
        calc_type += "_"

        # ensure separate folders are in place for each calculation input
        counter = 0
        # Count .traj outputs
        opt_restarts = 0

        # check/make folders
        while calc_type + filename + "_" + str(counter) \
                in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, calc_type + "*")]:
            counter += 1

        subdirectory_name = calc_type + filename + "_" + str(counter)
        subdirectory_name_prev = calc_type + filename + "_" + str(counter - 1)

        # check previous calculations for convergence
        if restart and counter > 0:
            if verbose:
                print("Previous calculation detected in", calc_type + filename + "_" + str(counter - 1))
            os.chdir(subdirectory_name_prev)

            if os.path.exists("evaluated_structures.traj"):
                self.prev_calcs = read("evaluated_structures.traj@:")
            elif verbose:
                print('The "evaluated_structures.traj" file not found, starting from scratch.')

            # check for number of restarted optimisations
            while str(counter - 1) + "_" + filename + "_" + str(opt_restarts) + ".traj" \
                    in [fn for fn in os.listdir() if fnmatch.fnmatch(fn, "*.traj")]:
                opt_restarts += 1

            # read the last optimisation
            if not internal:
                if os.path.exists(str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj"):
                    self.initial = read(str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj")
                    if verbose:
                        print("Restarting from", str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj")
            else:
                if os.path.exists(str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj"):
                    self.initial = read(str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj")
                    if verbose:
                        print("Restarting from", str(counter - 1) + "_" + filename + "_" + str(opt_restarts - 1) + ".traj")
                elif os.path.exists("geometry.in.next_step"):
                    shutil.copy("geometry.in.next_step", "next_step.in")
                    self.initial = read("next_step.in")
                    if verbose:
                        print("Restarting from geometry.in.next_step")

            os.chdir("..")

        return counter, subdirectory_name



    """
    def serialize(self):
        '''Save the instance to a file'''
        pass

    def recover(self):
        '''Recover a saved instance form a file'''
    """

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
    if sockets and not internal:
        from carmm.run.aims_calculator import get_aims_and_sockets_calculator
        # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
        # we need to specifically state what the name of the login node is so the two packages can communicate
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=dimensions,
                                                                 verbose=True,
                                                                 codata_warning=False
                                                                 )
    else:
        from carmm.run.aims_calculator import get_aims_calculator
        # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
        # we need to specifically state what the name of the login node is so the two packages can communicate
        fhi_calc = get_aims_calculator(dimensions=dimensions)

    # remove previous xc argument to ensure libxc warning override is first
    fhi_calc.parameters.pop("xc")
    fhi_calc.set(override_warning_libxc='True')

    # Forces required for optimisation:
    if not forces:
        fhi_calc.parameters.pop("compute_forces")

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

    if sockets and not internal:
        return sockets_calc, fhi_calc
    else:
        return fhi_calc


"""
    def search_ts_dyneb(self, initial, final, fmax, interpolation="idpp", restart=True, input_check=0.01, verbose=True):

        '''
        TODO: This requires parallelisation over images, otherwise with a shared calculator EVERY structure in the chain
            is calculated again per iteration. Waste of resources! This might be a bug, since the calculations were
            supposed to be dynamic. Need to investigate - maybe then the non-ML NEBS can be viable.
        Args:
            initial: 
            final: 
            fmax: 
            interpolation: 
            restart: 
            input_check: 
            verbose: 

        Returns:

        '''
        from ase.dyneb import DyNEB
        from ase.optimize import BFGS
        '''Retrieve common properties'''
        basis_set = self.basis_set
        dimensions = self.dimensions
        hpc = self.hpc
        params = self.params

        # Set the environment parameters
        set_aims_command(hpc=hpc, basis_set=basis_set)

        if not interpolation:
            interpolation = "idpp"

        # Read the geometry
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()

        # TODO: restart dedicated to dyneb
        counter, subdirectory_name = self._restart_setup(self.initial, filename, restart=restart,
                                                                   verbose=verbose)
        out = str(counter) + "_" + str(filename) + ".out"

        # ensure input is converged
        if input_check:  # TODO naming scheme + check for these files or initial/final.traj

            if not is_converged(initial, input_check):
                initial = self.aims_optimise(initial, input_check, restart=False, verbose=False)[0]
            if not is_converged(final, input_check):
                final = self.aims_optimise(final, input_check, restart=False, verbose=False)[0]

        if not os.path.exists(subdirectory_name):
            os.mkdir(subdirectory_name)
        os.chdir(subdirectory_name)

        # Make a band consisting of n images:
        n = 10 # TODO: do not hardcode
        images = [initial]
        images += [initial.copy() for i in range(n-1)]
        images += [final]
        neb = DyNEB(images,
                    k=0.5,
                    fmax=0.05,
                    climb=True,
                    parallel=False,
                    remove_rotation_and_translation=False,
                    world=None,
                    dynamic_relaxation=True,
                    scale_fmax=0.5,
                    method='aseneb',
                    allow_shared_calculator= True,
                    )
        # Interpolate linearly the positions of the three middle images:
        neb.interpolate(method=interpolation, mic=True, apply_constraint=True)

        # Create the sockets calculator - using a with statement means the object is closed at the end.
        with _calc_generator(params, out_fn=out, dimensions=dimensions)[0] as calculator:
            # Set calculators:
            for image in images[1:n]:
                image.calc = calculator

            optimizer = BFGS(neb,
                               trajectory=str(counter) + "_NEB_" + filename + ".traj",
                               )
            optimizer.run(fmax=0.05)

        # Find maximum energy, i.e. transition state to return it
        neb = read(str(counter) + "_NEB_" + filename + ".traj@-" +str(n) +":")
        self.ts = sorted(neb, key=lambda k: k.get_potential_energy(), reverse=True)[0]
        os.chdir("..")
"""

