# Author: Igor Kowalec
from ase import Atoms
from ase.io import read
from carmm.run.workflows.helper import CalculationHelper

class ReactMACE:
    """
    Class for streamlining calculations in an ASE/MACE setup.
    This workflow uses MACE-MP: Materials Project Force Fields

    Cite: 	arXiv:2401.00096 [physics.chem-ph]

    To install MACE via pip:
        > pip install --upgrade pip
        > pip install torch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 --index-url https://download.pytorch.org/whl/cu118
        > pip install mace-torch

    If you want to make use of dispersion correction install also:
        > pip install torch-dftd
    """

    def __init__(self,
                 params: dict,
                 force_field: str = "mace_mp",
                 filename: str = None,
                 dry_run: bool = False,
                 verbose: bool = True):
        """
        Args:
            params: dict
                Dictionary containing keywords and parameters for running the MACE calculator
            force_field: str
                Available currently are pre-trained "mace_mp", "mace_off", "mace_anicc" machine-learned force fields.
            filename: str
                Naming convention to follow for subsequent calculations and restarts.
                If None the chemical formula is used.
            dry_run: bool
                A dry run flag for CI-test
            verbose: bool
                Enables the verbosity of the performed operations

        Returns ReactMACE object
        """

        """Define basic parameters"""
        self.data = {}
        self.params = params
        self.force_field = force_field.lower()
        assert self.force_field in ["mace_mp", "mace_off", "mace_anicc"], \
            'Please use the supported MACE models -"mace_mp", "mace_off" or "mace_anicc".'
        self.filename = filename
        self.verbose = verbose

        """Define additional parameters"""
        self.initial = None  # input for optimisation or input for NEB initial image
        self.model_optimised = None  # optimised geometry with calculator attached
        self.model_post_processed = None  # post processed geometry with new calculator attached
        self.final = None  # input final image for NEB
        self.ts = None  # TS geometry from NEB
        self.prev_calcs = None  # for NEB restart
        self.interpolation = None  # For NEB interpolation
        self.dimensions = None  # for

        """ Set the test flag"""
        self.dry_run = dry_run

    def _get_mace_calculator(self):
        """
        This function returns one of the supported MACE calculators or the EMT calculator if the
        dry_run flag is enabled.
        """
        from ase.calculators.emt import EMT
        from mace.calculators import mace_mp, mace_off, mace_anicc

        force_fields = {"mace_mp": mace_mp,
                        "mace_off": mace_off,
                        "mace_anicc": mace_anicc}
        if not self.dry_run:
            return force_fields[self.force_field](**self.params)
        else:
            return EMT()

    def mace_optimise(self, atoms: Atoms, fmax=0.05, restart=True, relax_unit_cell=False):
        """
        Args:
            atoms: Atoms object
                Molecular structure to be optimised
            fmax: float
                Force convergence criterion in eV/Å
            restart: bool
                If True, seeks previously converged calculations based on filename
            relax_unit_cell: bool
                Request unit cell relaxation for periodic bulk structures

        Returns: Atoms object
            Optimised structure
        """

        import os

        self._initialize_parameters(atoms)
        helper = CalculationHelper(calc_type="Opt",
                                   parent_dir=os.getcwd(),
                                   filename=self.filename,
                                   restart=restart,
                                   verbose=self.verbose)

        counter, out, subdirectory_name, initial = helper.restart_setup()

        if initial:
            self.initial = initial

        self._perform_optimisation(subdirectory_name, counter, fmax, relax_unit_cell)

        return self.model_optimised

    def _perform_optimisation(self, subdirectory_name: str, counter: int, fmax: float,
                              relax_unit_cell: bool):
        """
        An internal function used in mace_optimise to resolve the working directory and perform the optimisation
        calculation of a given structure.

        Args:
            subdirectory_name: string
            counter: int
            fmax: float
            relax_unit_cell: bool

        Returns:None

        """
        from ase.optimize import BFGS
        from carmm.analyse.forces import is_converged
        import os


        opt_restarts = 0

        if not is_converged(self.initial, fmax):
            os.makedirs(subdirectory_name, exist_ok=True)
            self.initial.calc = self._get_mace_calculator()

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

    def search_ts_neb(self, initial, final, fmax, n, k=0.05, method="aseneb", interpolation="idpp", input_check=0.01,
                      max_steps=100, restart=True):
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
            k: float or list of floats
                Spring constant(s) in eV/Ang.  One number or one for each spring.
            method: str
                NEB method for the CI-NEB as implemented in ASE, 'string' by default
            interpolation: str or []
                The "idpp" or "linear" interpolation types are supported in ASE. alternatively user can provide a
                custom interpolation as a list of Atoms objects.
            input_check: float or None
                If float the calculators of the input structures will be checked if the structures are below
                the requested fmax and an optimisation will be performed if not.
            max_steps: int
                Maximum number of iteration before stopping the optimizer
            restart: bool
                If True, seeks previously converged calculations based on filename

        Returns: Atoms object
            Transition state geometry structure
        """

        from ase.neb import NEB
        from ase.optimize import FIRE
        from carmm.analyse.forces import is_converged
        import os

        '''Retrieve common properties'''
        parent_dir = os.getcwd()

        '''Read the geometry and establish a naming convention'''
        if self.filename:
            filename = self.filename
        else:
            filename = initial.get_chemical_formula()

        '''Ensure input is converged'''
        if input_check:
            if not is_converged(initial, input_check):
                self.filename = filename + "_initial"
                initial = self.mace_optimise(initial, fmax=input_check, restart=True)
            if not is_converged(final, input_check):
                self.filename = filename + "_final"
                final = self.mace_optimise(final, fmax=input_check, restart=True)

            '''Revert to original name'''
            self.filename = filename

        '''Setup the TS calculation'''
        helper = CalculationHelper(calc_type="TS", parent_dir=os.getcwd(), filename=self.filename, restart=restart,
                                   verbose=self.verbose)
        counter, out, subdirectory_name, minimum_energy_path = helper.restart_setup()

        if not minimum_energy_path:
            minimum_energy_path = [None, None]
        elif restart:
            mep = minimum_energy_path[1][-n:]
            previous_neb = NEB(mep, k=k, method=method, climb=True, parallel=True, allow_shared_calculator=False)

            previous_neb.get_forces()
            residual = previous_neb.get_residual()

            if residual <= fmax:
                self.interpolation = mep
                minimum_energy_path[0] = mep
            else:
                self.interpolation = mep
                interpolation = mep

        if not minimum_energy_path[0]:
            """Create the calculators for all images"""

            os.makedirs(subdirectory_name, exist_ok=True)
            os.chdir(subdirectory_name)

            if interpolation in ["idpp", "linear"]:
                images = [initial]
                for i in range(n):
                    image = initial.copy()
                    image.calc = self._get_mace_calculator()
                    images.append(image)
                images.append(final)

            elif isinstance(interpolation, list):
                assert [isinstance(i, Atoms) for i in interpolation], \
                    "Interpolation must be a list of Atoms objects, 'idpp' or 'linear'!"
                images = interpolation
                for i in range(1, len(interpolation) - 1):
                    images[i].calc = self._get_mace_calculator()
            else:
                raise ValueError("Interpolation must be a list of Atoms objects, 'idpp' or 'linear'!")

            neb = NEB(images, k=k, method=method, climb=True, parallel=True, allow_shared_calculator=False)

            if interpolation in ["idpp", "linear"]:
                neb.interpolate(method=interpolation, mic=True, apply_constraint=True)

            qn = FIRE(neb, trajectory=f'{self.filename}_NEB.traj')
            qn.run(fmax=fmax, steps=max_steps)

            minimum_energy_path[0] = images
            self.interpolation = minimum_energy_path[0]

        '''Find maximum energy, i.e. transition state to return it'''
        self.ts = sorted(minimum_energy_path[0], key=lambda x: x.get_potential_energy(), reverse=True)[0]
        os.chdir(parent_dir)

        return self.ts
