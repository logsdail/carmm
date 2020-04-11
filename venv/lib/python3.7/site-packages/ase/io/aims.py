import time
import warnings

from ase.units import Ang, fs

v_unit = Ang / (1000.0 * fs)


def read_aims(filename, apply_constraints=True):
    """Import FHI-aims geometry type files.

    Reads unitcell, atom positions and constraints from
    a geometry.in file.

    If geometric constraint (symmetry parameters) are in the file
    include that information in atoms.info["symmetry_block"]
    """

    from ase import Atoms
    from ase.constraints import (
        FixAtoms,
        FixCartesian,
        FixScaledParametricRelations,
        FixCartesianParametricRelations,
    )
    import numpy as np

    atoms = Atoms()
    with open(filename, "r") as fd:
        lines = fd.readlines()

    positions = []
    cell = []
    symbols = []
    velocities = []
    magmoms = []
    symmetry_block = []
    charges = []
    fix = []
    fix_cart = []
    xyz = np.array([0, 0, 0])
    i = -1
    n_periodic = -1
    periodic = np.array([False, False, False])
    cart_positions, scaled_positions = False, False
    for line in lines:
        inp = line.split()
        if inp == []:
            continue
        if inp[0] == "atom":
            cart_positions = True
            if xyz.all():
                fix.append(i)
            elif xyz.any():
                fix_cart.append(FixCartesian(i, xyz))
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            positions.append(floatvect)
            magmoms.append(0.0)
            charges.append(0.0)
            symbols.append(inp[-1])
            i += 1
            xyz = np.array([0, 0, 0])
        elif inp[0] == "atom_frac":
            scaled_positions = True
            if xyz.all():
                fix.append(i)
            elif xyz.any():
                fix_cart.append(FixCartesian(i, xyz))
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            positions.append(floatvect)
            magmoms.append(0.)
            symbols.append(inp[-1])
            i += 1
            xyz = np.array([0, 0, 0])

        elif inp[0] == "lattice_vector":
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            cell.append(floatvect)
            n_periodic = n_periodic + 1
            periodic[n_periodic] = True

        elif inp[0] == "initial_moment":
            magmoms[-1] = float(inp[1])

        elif inp[0] == "initial_charge":
            charges[-1] = float(inp[1])

        elif inp[0] == "constrain_relaxation":
            if inp[1] == ".true.":
                fix.append(i)
            elif inp[1] == "x":
                xyz[0] = 1
            elif inp[1] == "y":
                xyz[1] = 1
            elif inp[1] == "z":
                xyz[2] = 1

        elif inp[0] == "velocity":
            floatvect = [v_unit * float(l) for l in inp[1:4]]
            velocities.append(floatvect)

        elif inp[0] in [
            "symmetry_n_params",
            "symmetry_params",
            "symmetry_lv",
            "symmetry_frac",
        ]:
            symmetry_block.append(" ".join(inp))

    if xyz.all():
        fix.append(i)
    elif xyz.any():
        fix_cart.append(FixCartesian(i, xyz))

    if cart_positions and scaled_positions:
        raise Exception(
            "Can't specify atom positions with mixture of "
            "Cartesian and fractional coordinates"
        )
    elif scaled_positions and periodic.any():
        atoms = Atoms(
            symbols, scaled_positions=positions, cell=cell, pbc=periodic
        )
    else:
        atoms = Atoms(symbols, positions)

    if len(velocities) > 0:
        if len(velocities) != len(positions):
            raise Exception(
                "Number of positions and velocities have to coincide."
            )
        atoms.set_velocities(velocities)

    fix_params = []

    if len(symmetry_block) > 5:
        params = symmetry_block[1].split()[1:]

        lattice_expressions = []
        lattice_params = []

        atomic_expressions = []
        atomic_params = []

        n_lat_param = int(symmetry_block[0].split(" ")[2])

        lattice_params = params[:n_lat_param]
        atomic_params = params[n_lat_param:]

        for ll, line in enumerate(symmetry_block[2:]):
            expression = " ".join(line.split(" ")[1:])
            if ll < 3:
                lattice_expressions += expression.split(",")
            else:
                atomic_expressions += expression.split(",")

        fix_params.append(
            FixCartesianParametricRelations.from_expressions(
                list(range(3)),
                lattice_params,
                lattice_expressions,
                use_cell=True,
            )
        )

        fix_params.append(
            FixScaledParametricRelations.from_expressions(
                list(range(len(atoms))), atomic_params, atomic_expressions
            )
        )

    if any(magmoms):
        atoms.set_initial_magnetic_moments(magmoms)
    if any(charges):
        atoms.set_initial_charges(charges)

    if periodic.any():
        atoms.set_cell(cell)
        atoms.set_pbc(periodic)
    if len(fix):
        atoms.set_constraint([FixAtoms(indices=fix)] + fix_cart + fix_params)
    else:
        atoms.set_constraint(fix_cart + fix_params)

    if fix_params and apply_constraints:
        atoms.set_positions(atoms.get_positions())
    return atoms


def write_aims(
    filename,
    atoms,
    scaled=False,
    geo_constrain=False,
    velocities=False,
    ghosts=None,
    info_str=None,
    wrap=False,
):
    """Method to write FHI-aims geometry files.

    Writes the atoms positions and constraints (only FixAtoms is
    supported at the moment).

    Args:
        filename: str
            Name of file to output structure to
        atoms: ase.atoms.Atoms
            structure to output to the file
        scaled: bool
            If True use fractional coordinates instead of Cartesian coordinates
        symmetry_block: list of str
            List of geometric constraints as defined in:
            https://arxiv.org/abs/1908.01610
        velocities: bool
            If True add the atomic velocity vectors to the file
        ghosts: list of Atoms
            A list of ghost atoms for the system
        info_str: str
            A string to be added to the header of the file
        wrap: bool
            Wrap atom positions to cell before writing
    """

    from ase.constraints import FixAtoms, FixCartesian

    import numpy as np

    if isinstance(atoms, (list, tuple)):
        if len(atoms) > 1:
            raise RuntimeError(
                "Don't know how to save more than "
                "one image to FHI-aims input"
            )
        else:
            atoms = atoms[0]

    if geo_constrain:
        if not scaled:
            warnings.warn(
                "Setting scaled to True because a symmetry_block is detected."
            )
            scaled = True

    fd = open(filename, "w")
    fd.write("#=======================================================\n")
    fd.write("# FHI-aims file: " + filename + "\n")
    fd.write("# Created using the Atomic Simulation Environment (ASE)\n")
    fd.write("# " + time.asctime() + "\n")

    # If writing additional information is requested via info_str:
    if info_str is not None:
        fd.write("\n# Additional information:\n")
        if isinstance(info_str, list):
            fd.write("\n".join(["#  {}".format(s) for s in info_str]))
        else:
            fd.write("# {}".format(info_str))
        fd.write("\n")

    fd.write("#=======================================================\n")

    i = 0
    if atoms.get_pbc().any():
        for n, vector in enumerate(atoms.get_cell()):
            fd.write("lattice_vector ")
            for i in range(3):
                fd.write("%16.16f " % vector[i])
            fd.write("\n")
    fix_cart = np.zeros([len(atoms), 3])

    # else aims crashes anyways
    # better be more explicit
    # write_magmoms = np.any([a.magmom for a in atoms])

    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixAtoms):
                fix_cart[constr.index] = [1, 1, 1]
            elif isinstance(constr, FixCartesian):
                fix_cart[constr.a] = -constr.mask + 1

    if ghosts is None:
        ghosts = np.zeros(len(atoms))
    else:
        assert len(ghosts) == len(atoms)

    if geo_constrain:
        wrap = False
    scaled_positions = atoms.get_scaled_positions(wrap=wrap)

    for i, atom in enumerate(atoms):
        if ghosts[i] == 1:
            atomstring = "empty "
        elif scaled:
            atomstring = "atom_frac "
        else:
            atomstring = "atom "
        fd.write(atomstring)
        if scaled:
            for pos in scaled_positions[i]:
                fd.write("%16.16f " % pos)
        else:
            for pos in atom.position:
                fd.write("%16.16f " % pos)
        fd.write(atom.symbol)
        fd.write("\n")
        # (1) all coords are constrained:
        if fix_cart[i].all():
            fd.write("    constrain_relaxation .true.\n")
        # (2) some coords are constrained:
        elif fix_cart[i].any():
            xyz = fix_cart[i]
            for n in range(3):
                if xyz[n]:
                    fd.write("    constrain_relaxation %s\n" % "xyz"[n])
        if atom.charge:
            fd.write("    initial_charge %16.6f\n" % atom.charge)
        if atom.magmom:
            fd.write("    initial_moment %16.6f\n" % atom.magmom)

        # Write velocities if this is wanted
        if velocities and atoms.get_velocities() is not None:
            fd.write(
                "    velocity {:.16f} {:.16f} {:.16f}\n".format(
                    *atoms.get_velocities()[i] / v_unit
                )
            )

    if geo_constrain:
        for line in get_sym_block(atoms):
            fd.write(line)


def get_sym_block(atoms):
    """Get the symmetry block for the Parametric constraints in atoms.constraints"""
    import numpy as np

    from ase.constraints import (
        FixScaledParametricRelations,
        FixCartesianParametricRelations,
    )

    # Initialize param/expressions lists
    atomic_sym_params = []
    lv_sym_params = []
    atomic_param_constr = np.zeros((len(atoms),), dtype="<U100")
    lv_param_constr = np.zeros((3,), dtype="<U100")

    # Populate param/expressions list
    for constr in atoms.constraints:
        if isinstance(constr, FixScaledParametricRelations):
            atomic_sym_params += constr.params

            if np.any(atomic_param_constr[constr.indices] != ""):
                warnings.warn(
                    "multiple parametric constraints defined for the same atom, using the last one defined"
                )

            atomic_param_constr[constr.indices] = [
                ", ".join(expression) for expression in constr.expressions
            ]
        elif isinstance(constr, FixCartesianParametricRelations):
            lv_sym_params += constr.params

            if np.any(lv_param_constr[constr.indices] != ""):
                warnings.warn(
                    "multiple parametric constraints defined for the same lattice vector, using the last one defined"
                )

            lv_param_constr[constr.indices] = [
                ", ".join(expression) for expression in constr.expressions
            ]

    if np.all(atomic_param_constr == "") and np.all(lv_param_constr == ""):
        return []

    # Check Constraint Parameters
    if len(atomic_sym_params) != len(np.unique(atomic_sym_params)):
        warnings.warn(
            "Some parameters were used across constraints, they will be combined in the aims calculations"
        )
        atomic_sym_params = np.unique(atomic_sym_params)

    if len(lv_sym_params) != len(np.unique(lv_sym_params)):
        warnings.warn(
            "Some parameters were used across constraints, they will be combined in the aims calculations"
        )
        lv_sym_params = np.unique(lv_sym_params)

    if np.any(atomic_param_constr == ""):
        raise IOError(
            "FHI-aims input files require all atoms have defined parametric constraints"
        )

    cell_inds = np.where(lv_param_constr == "")[0]
    for ind in cell_inds:
        lv_param_constr[ind] = "{:.16f}, {:.16f}, {:.16f}".format(
            *atoms.cell[ind]
        )

    n_atomic_params = len(atomic_sym_params)
    n_lv_params = len(lv_sym_params)
    n_total_params = n_atomic_params + n_lv_params

    sym_block = []
    if n_total_params > 0:
        sym_block.append(
            "#=======================================================\n"
        )
        sym_block.append("# Parametric constraints\n")
        sym_block.append(
            "#=======================================================\n"
        )
        sym_block.append(
            "symmetry_n_params {:d} {:d} {:d}\n".format(
                n_total_params, n_lv_params, n_atomic_params
            )
        )
        sym_block.append(
            "symmetry_params %s\n"
            % " ".join(lv_sym_params + atomic_sym_params)
        )

        for constr in lv_param_constr:
            sym_block.append("symmetry_lv {:s}\n".format(constr))

        for constr in atomic_param_constr:
            sym_block.append("symmetry_frac {:s}\n".format(constr))
    return sym_block


# except KeyError:
#     continue


def read_energy(filename):
    for line in open(filename, "r"):
        if line.startswith("  | Total energy corrected"):
            E = float(line.split()[-2])
    return E


def read_aims_output(filename, index=-1):
    """Import FHI-aims output files with all data available, i.e.
    relaxations, MD information, force information etc etc etc."""
    from ase import Atoms, Atom
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.constraints import FixAtoms, FixCartesian

    molecular_dynamics = False
    fd = open(filename, "r")
    cell = []
    images = []
    fix = []
    fix_cart = []
    f = None
    pbc = False
    found_aims_calculator = False
    stress = None
    for line in fd:
        # if "List of parameters used to initialize the calculator:" in line:
        #     next(fd)
        #     calc = read_aims_calculator(fd)
        #     calc.out = filename
        #     found_aims_calculator = True
        if "| Number of atoms                   :" in line:
            inp = line.split()
            n_atoms = int(inp[5])
        if "| Unit cell:" in line:
            if not pbc:
                pbc = True
                for i in range(3):
                    inp = next(fd).split()
                    cell.append([inp[1], inp[2], inp[3]])
        if "Found relaxation constraint for atom" in line:
            xyz = [0, 0, 0]
            ind = int(line.split()[5][:-1]) - 1
            if "All coordinates fixed" in line:
                if ind not in fix:
                    fix.append(ind)
            if "coordinate fixed" in line:
                coord = line.split()[6]
                if coord == "x":
                    xyz[0] = 1
                elif coord == "y":
                    xyz[1] = 1
                elif coord == "z":
                    xyz[2] = 1
                keep = True
                for n, c in enumerate(fix_cart):
                    if ind == c.a:
                        keep = False
                if keep:
                    fix_cart.append(FixCartesian(ind, xyz))
                else:
                    fix_cart[n].mask[xyz.index(1)] = 0
        if "Atomic structure:" in line and not molecular_dynamics:
            next(fd)
            atoms = Atoms()
            for i in range(n_atoms):
                inp = next(fd).split()
                atoms.append(Atom(inp[3], (inp[4], inp[5], inp[6])))
        if "Complete information for previous time-step:" in line:
            molecular_dynamics = True
        if "Updated atomic structure:" in line and not molecular_dynamics:
            next(fd)
            atoms = Atoms()
            for i in range(n_atoms):
                inp = next(fd).split()
                if "lattice_vector" in inp[0]:
                    cell = []
                    for i in range(3):
                        cell += [[float(inp[1]), float(inp[2]), float(inp[3])]]
                        inp = next(fd).split()
                    atoms.set_cell(cell)
                    inp = next(fd).split()
                atoms.append(Atom(inp[4], (inp[1], inp[2], inp[3])))
                if molecular_dynamics:
                    inp = next(fd).split()
        if "Atomic structure (and velocities)" in line:
            next(fd)
            atoms = Atoms()
            velocities = []
            for i in range(n_atoms):
                inp = next(fd).split()
                atoms.append(Atom(inp[4], (inp[1], inp[2], inp[3])))
                inp = next(fd).split()
                floatvect = [v_unit * float(l) for l in inp[1:4]]
                velocities.append(floatvect)
            atoms.set_velocities(velocities)
            if len(fix):
                atoms.set_constraint([FixAtoms(indices=fix)] + fix_cart)
            else:
                atoms.set_constraint(fix_cart)
            images.append(atoms)

        # FlK: add analytical stress and replace stress=None
        if "Analytical stress tensor - Symmetrized" in line:
            # scroll to significant lines
            for _ in range(4):
                next(fd)
            stress = []
            for _ in range(3):
                inp = next(fd)
                stress.append([float(i) for i in inp.split()[2:5]])

        if "Total atomic forces" in line:
            f = []
            for i in range(n_atoms):
                inp = next(fd).split()
                # FlK: use inp[-3:] instead of inp[1:4] to make sure this works
                # when atom number is not preceded by a space.
                f.append([float(i) for i in inp[-3:]])
            if not found_aims_calculator:
                e = images[-1].get_potential_energy()
                # FlK: Add the stress if it has been computed
                if stress is None:
                    calc = SinglePointCalculator(atoms, energy=e, forces=f)
                else:
                    calc = SinglePointCalculator(
                        atoms, energy=e, forces=f, stress=stress
                    )
                images[-1].set_calculator(calc)
            e = None
            f = None

        if "Total energy corrected" in line:
            e = float(line.split()[5])
            if pbc:
                atoms.set_cell(cell)
                atoms.pbc = True
            if not found_aims_calculator:
                atoms.set_calculator(SinglePointCalculator(atoms, energy=e))
            if not molecular_dynamics:
                if len(fix):
                    atoms.set_constraint([FixAtoms(indices=fix)] + fix_cart)
                else:
                    atoms.set_constraint(fix_cart)
                images.append(atoms)
            e = None
            # if found_aims_calculator:
            # calc.set_results(images[-1])
            # images[-1].set_calculator(calc)

        # FlK: add stress per atom
        if "Per atom stress (eV) used for heat flux calculation" in line:
            # scroll to boundary
            next(l for l in fd if "-------------" in l)

            stresses = []
            for l in [next(fd) for _ in range(n_atoms)]:
                # Read stresses
                xx, yy, zz, xy, xz, yz = [float(d) for d in l.split()[2:8]]
                stresses.append([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])

            if not found_aims_calculator:
                e = images[-1].get_potential_energy()
                f = images[-1].get_forces()
                stress = images[-1].get_stress(voigt=False)

                calc = SinglePointCalculator(
                    atoms, energy=e, forces=f, stress=stress, stresses=stresses
                )
                images[-1].set_calculator(calc)

    fd.close()
    if molecular_dynamics:
        images = images[1:]

    # return requested images, code borrowed from ase/io/trajectory.py
    if isinstance(index, int):
        return images[index]
    else:
        step = index.step or 1
        if step > 0:
            start = index.start or 0
            if start < 0:
                start += len(images)
            stop = index.stop or len(images)
            if stop < 0:
                stop += len(images)
        else:
            if index.start is None:
                start = len(images) - 1
            else:
                start = index.start
                if start < 0:
                    start += len(images)
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop
                if stop < 0:
                    stop += len(images)
        return [images[i] for i in range(start, stop, step)]
