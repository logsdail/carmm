"""
ASE Calculator for interatomic models compatible with the Knowledgebase
of Interatomic Models (KIM) application programming interface (API).
Written by:

Mingjian Wen
Daniel S. Karls
University of Minnesota
"""
import numpy as np

from ase.calculators.calculator import Calculator
from ase.calculators.calculator import compare_atoms

from . import kimpy_wrappers
from . import neighborlist


class KIMModelData(object):
    """Initializes and subsequently stores the KIM API Portable Model
    object, KIM API ComputeArguments object, and the neighbor list
    object used by instances of KIMModelCalculator.  Also stores the
    arrays which are registered in the KIM API and which are used to
    communicate with the model.
    """

    def __init__(self, model_name, ase_neigh, neigh_skin_ratio, debug=False):
        self.model_name = model_name
        self.ase_neigh = ase_neigh
        self.debug = debug

        # Initialize KIM API Portable Model object and ComputeArguments object
        self.init_kim()

        # Set cutoff
        model_influence_dist = self.kim_model.get_influence_distance()
        model_cutoffs, padding_not_require_neigh = (
            self.kim_model.get_neighbor_list_cutoffs_and_hints()
        )

        self.species_map = self.create_species_map()

        # Initialize neighbor list object
        self.init_neigh(
            neigh_skin_ratio,
            model_influence_dist,
            model_cutoffs,
            padding_not_require_neigh,
        )

    def __del__(self):
        self.clean()

    def init_kim(self):
        """Create the KIM API Portable Model object and KIM API ComputeArguments
        object
        """
        if self.kim_initialized:
            return

        self.kim_model = kimpy_wrappers.PortableModel(self.model_name, self.debug)

        # KIM API model object is what actually creates/destroys the ComputeArguments
        # object, so we must pass it as a parameter
        self.compute_args = self.kim_model.compute_arguments_create()

    def init_neigh(
        self,
        neigh_skin_ratio,
        model_influence_dist,
        model_cutoffs,
        padding_not_require_neigh,
    ):
        """Initialize neighbor list, either an ASE-native neighborlist
        or one created using the neighlist module in kimpy
        """
        neigh_list_object_type = (
            neighborlist.ASENeighborList
            if self.ase_neigh
            else neighborlist.KimpyNeighborList
        )
        self.neigh = neigh_list_object_type(
            self.compute_args,
            neigh_skin_ratio,
            model_influence_dist,
            model_cutoffs,
            padding_not_require_neigh,
            self.debug,
        )

    def update_compute_args_pointers(self, energy, forces):
        self.compute_args.update(
            self.num_particles,
            self.species_code,
            self.particle_contributing,
            self.coords,
            energy,
            forces,
        )

    def create_species_map(self):
        """Get all the supported species of the KIM model and the
        corresponding integer codes used by the model

        Returns
        -------
        species_map : dict
            key : str
                chemical symbols (e.g. "Ar")
            value : int
                species integer code (e.g. 1)
        """
        supported_species, codes = self.get_model_supported_species_and_codes()
        species_map = dict()
        for i, spec in enumerate(supported_species):
            species_map[spec] = codes[i]
            if self.debug:
                print(
                    "Species {} is supported and its code is: {}".format(spec, codes[i])
                )

        return species_map

    def clean_neigh(self):
        """If the neighbor list method being used is the one in the
        kimpy neighlist module, deallocate its memory
        """
        if self.neigh_initialized:
            self.neigh.clean()
            del self.neigh

    def clean_kim(self):
        """Deallocate the memory allocated to the KIM API Portable Model object
        and KIM API ComputeArguments object
        """
        if self.kim_initialized:
            self.kim_model.compute_arguments_destroy(self.compute_args)
            self.kim_model.destroy()
            del self.kim_model

    def clean(self):
        """Deallocate the KIM API Portable Model object, KIM API ComputeArguments
        object, and, if applicable, the neighbor list object
        """
        self.clean_neigh()
        self.clean_kim()

    @property
    def padding_image_of(self):
        return self.neigh.padding_image_of

    @property
    def num_particles(self):
        return self.neigh.num_particles

    @property
    def coords(self):
        return self.neigh.coords

    @property
    def particle_contributing(self):
        return self.neigh.particle_contributing

    @property
    def species_code(self):
        return self.neigh.species_code

    @property
    def kim_initialized(self):
        return hasattr(self, "kim_model")

    @property
    def neigh_initialized(self):
        return hasattr(self, "neigh")

    @property
    def get_model_supported_species_and_codes(self):
        return self.kim_model.get_model_supported_species_and_codes


class KIMModelCalculator(Calculator):
    """Calculator that works with KIM Portable Models (PMs).

    Calculator that carries out direct communication between ASE and a
    KIM Portable Model (PM) through the kimpy library (which provides a
    set of python bindings to the KIM API).

    Parameters
    ----------
    model_name : str
      The unique identifier assigned to the interatomic model (for
      details, see https://openkim.org/doc/schema/kim-ids)

    ase_neigh : bool, optional
      False (default): Use kimpy's neighbor list library

      True: Use ASE's internal neighbor list mechanism (usually slower
      than the kimpy neighlist library)

    neigh_skin_ratio : float, optional
      Used to determine the neighbor list cutoff distance, r_neigh,
      through the relation r_neigh = (1 + neigh_skin_ratio) * rcut,
      where rcut is the model's influence distance. (Default: 0.2)

    release_GIL : bool, optional
      Whether to release python GIL.  Releasing the GIL allows a KIM
      model to run with multiple concurrent threads. (Default: False)

    debug : bool, optional
      If True, detailed information is printed to stdout. (Default:
      False)
    """

    implemented_properties = ["energy", "forces", "stress"]

    def __init__(
        self,
        model_name,
        ase_neigh=False,
        neigh_skin_ratio=0.2,
        release_GIL=False,
        debug=False,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)

        self.model_name = model_name
        self.release_GIL = release_GIL
        self.debug = debug

        if neigh_skin_ratio < 0:
            raise ValueError('Argument "neigh_skin_ratio" must be non-negative')

        # Model output
        self.energy = None
        self.forces = None

        # Create KIMModelData object. This will take care of creating and storing the KIM
        # API Portable Model object, KIM API ComputeArguments object, and the neighbor
        # list object that our calculator needs
        self.kimmodeldata = KIMModelData(
            self.model_name, ase_neigh, neigh_skin_ratio, self.debug
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, value, traceback):
        # Explicitly deallocate all three objects held by the KIMModelData
        # instance referenced by our calculator
        self.kimmodeldata.clean()

    def __repr__(self):
        return "KIMModelCalculator(model_name={})".format(self.model_name)

    def calculate(
        self,
        atoms=None,
        properties=["energy", "forces", "stress"],
        system_changes=["positions", "numbers", "cell", "pbc"],
    ):
        """
        Inherited method from the ase Calculator class that is called by
        get_property()

        Parameters
        ----------
        atoms : Atoms
            Atoms object whose properties are desired

        properties : list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces' and 'stress'.

        system_changes : list of str
            List of what has changed since last calculation.  Can be any
            combination of these six: 'positions', 'numbers', 'cell',
            and 'pbc'.
        """

        Calculator.calculate(self, atoms, properties, system_changes)

        # Update KIM API input data and neighbor list, if necessary
        if system_changes:
            if self.need_neigh_update(atoms, system_changes):
                self.update_neigh(atoms, self.species_map)
                self.energy = np.array([0.0], dtype=np.double)
                self.forces = np.zeros([self.num_particles[0], 3], dtype=np.double)
                self.update_compute_args_pointers(self.energy, self.forces)
            else:
                self.update_kim_coords(atoms)

            self.kim_model.compute(self.compute_args, self.release_GIL)

        energy = self.energy[0]
        forces = self.assemble_padding_forces()

        try:
            volume = atoms.get_volume()
            stress = self.compute_virial_stress(self.forces, self.coords, volume)
        except ValueError:  # Volume cannot be computed
            stress = None

        # Quantities passed back to ASE
        self.results["energy"] = energy
        self.results["free_energy"] = energy
        self.results["forces"] = forces
        self.results["stress"] = stress

    def check_state(self, atoms, tol=1e-15):
        return compare_atoms(self.atoms, atoms, excluded_properties={'initial_charges',
            'initial_magmoms'})

    def assemble_padding_forces(self):
        """
        Assemble forces on padding atoms back to contributing atoms.

        Parameters
        ----------
        forces : 2D array of doubles
            Forces on both contributing and padding atoms

        num_contrib:  int
            Number of contributing atoms

        padding_image_of : 1D array of int
            Atom number, of which the padding atom is an image


        Returns
        -------
            Total forces on contributing atoms.
        """

        total_forces = np.array(self.forces[: self.num_contributing_particles])

        if self.padding_image_of.size != 0:
            pad_forces = self.forces[self.num_contributing_particles :]
            for f, org_index in zip(pad_forces, self.padding_image_of):
                total_forces[org_index] += f

        return total_forces

    @staticmethod
    def compute_virial_stress(forces, coords, volume):
        """Compute the virial stress in Voigt notation.

        Parameters
        ----------
        forces : 2D array
            Partial forces on all atoms (padding included)

        coords : 2D array
            Coordinates of all atoms (padding included)

        volume : float
            Volume of cell

        Returns
        -------
        stress : 1D array
            stress in Voigt order (xx, yy, zz, yz, xz, xy)
        """
        stress = np.zeros(6)
        stress[0] = -np.dot(forces[:, 0], coords[:, 0]) / volume
        stress[1] = -np.dot(forces[:, 1], coords[:, 1]) / volume
        stress[2] = -np.dot(forces[:, 2], coords[:, 2]) / volume
        stress[3] = -np.dot(forces[:, 1], coords[:, 2]) / volume
        stress[4] = -np.dot(forces[:, 0], coords[:, 2]) / volume
        stress[5] = -np.dot(forces[:, 0], coords[:, 1]) / volume

        return stress

    def get_model_supported_species_and_codes(self):
        return self.kimmodeldata.get_model_supported_species_and_codes

    @property
    def update_compute_args_pointers(self):
        return self.kimmodeldata.update_compute_args_pointers

    @property
    def kim_model(self):
        return self.kimmodeldata.kim_model

    @property
    def compute_args(self):
        return self.kimmodeldata.compute_args

    @property
    def num_particles(self):
        return self.kimmodeldata.num_particles

    @property
    def coords(self):
        return self.kimmodeldata.coords

    @property
    def padding_image_of(self):
        return self.kimmodeldata.padding_image_of

    @property
    def species_map(self):
        return self.kimmodeldata.species_map

    @property
    def neigh(self):
        return self.kimmodeldata.neigh

    @property
    def num_contributing_particles(self):
        return self.neigh.num_contributing_particles

    @property
    def update_kim_coords(self):
        return self.neigh.update_kim_coords

    @property
    def need_neigh_update(self):
        return self.neigh.need_neigh_update

    @property
    def update_neigh(self):
        return self.neigh.update
