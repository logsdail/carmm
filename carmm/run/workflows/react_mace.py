# Author: Igor Kowalec
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read
from ase.optimize import BFGS
from carmm.run.workflows.helper import CalculationHelper
from carmm.analyse.forces import is_converged
from mace.calculators import mace_mp as calculator
import os

class ReactMACE:
    '''
    Class for streamlining calculations in an ASE/MACE setup.
    This workflow uses MACE-MP: Materials Project Force Fields

    Cite: 	arXiv:2401.00096 [physics.chem-ph]

    To install MACE via pip:
        > pip install --upgrade pip
        > pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
        > pip install mace-torch

    If you want to make use of dispersion correction install also:
        > pip install torch-dftd

    When installed, uncomment the import in line 8
    '''

    def __init__(self,
                 params: dict,
                 filename: str = None,
                 dry_run: bool = False,
                 verbose: bool = True):
        """
        Args:
            params: dict
                Dictionary containing keywords and parameters for running the MACE calculator
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
        self.filename = filename
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

    def mace_optimise(self, atoms: Atoms, fmax=0.05, restart=True, relax_unit_cell=False):

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

        return self.model_optimised


    def _perform_optimization(self, subdirectory_name: str, out: str, counter: int, fmax: float,
                              relax_unit_cell: bool):
        """
        An internal function used in mace_optimise to resolve the working directory and perform the optimisation
        calculation of a given structure.

        Args:
            subdirectory_name: string
            out: str
            counter: int
            fmax: float
            relax_unit_cell: bool

        Returns:None

        """
        global calculator
        opt_restarts = 0

        if not is_converged(self.initial, fmax):
            os.makedirs(subdirectory_name, exist_ok=True)

            if self.dry_run:
                calculator = EMT

            if not self.dry_run:
                self.initial.calc = calculator(**self.params)
            else:
                self.initial.calc = calculator()

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