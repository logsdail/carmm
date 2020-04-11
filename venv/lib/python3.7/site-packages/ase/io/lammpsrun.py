import gzip
import struct
from os.path import splitext
from collections import deque
import numpy as np

from ase.atoms import Atoms
from ase.quaternions import Quaternions
from ase.calculators.singlepoint import SinglePointCalculator
from ase.parallel import paropen
from ase.utils import basestring
from ase.calculators.lammps import convert


def read_lammps_dump(infileobj, **kwargs):
    """Method which reads a LAMMPS dump file.

       LAMMPS chooses output method depending on the given suffix:
        - .bin  : binary file
        - .gz   : output piped through gzip
        - .mpiio: using mpiio (should be like cleartext,
                  with different ordering)
        - else  : normal clear-text format

    :param infileobj: string to file, opened file or file-like stream

    """
    # !TODO: add support for lammps-regex naming schemes (output per
    # processor and timestep wildcards)

    if isinstance(infileobj, basestring):
        suffix = splitext(infileobj)[-1]
        if suffix == ".bin":
            fileobj = paropen(infileobj, "rb")
        elif suffix == ".gz":
            # !TODO: save for parallel execution?
            fileobj = gzip.open(infileobj, "rb")
        else:
            fileobj = paropen(infileobj)
    else:
        suffix = splitext(infileobj.name)[-1]
        fileobj = infileobj

    if suffix == ".bin":
        return read_lammps_dump_binary(fileobj, **kwargs)

    return read_lammps_dump_text(fileobj, **kwargs)


def lammps_data_to_ase_atoms(
    data,
    colnames,
    cell,
    celldisp,
    pbc=False,
    atomsobj=Atoms,
    order=True,
    specorder=None,
    prismobj=None,
    units="metal",
):
    """Extract positions and other per-atom parameters and create Atoms

    :param data: per atom data
    :param colnames: index for data
    :param cell: cell dimensions
    :param celldisp: origin shift
    :param pbc: periodic boundaries
    :param atomsobj: function to create ase-Atoms object
    :param order: sort atoms by id. Might be faster to turn off
    :param specorder: list of species to map lammps types to ase-species
    (usually .dump files to not contain type to species mapping)
    :param prismobj: Coordinate transformation between lammps and ase
    :type prismobj: Prism
    :param units: lammps units for unit transformation between lammps and ase
    :returns: Atoms object
    :rtype: Atoms

    """
    # data array of doubles
    ids = data[:, colnames.index("id")].astype(int)
    types = data[:, colnames.index("type")].astype(int)
    if order:
        sort_order = np.argsort(ids)
        ids = ids[sort_order]
        data = data[sort_order, :]
        types = types[sort_order]

    # reconstruct types from given specorder
    if specorder:
        types = [specorder[t - 1] for t in types]

    def get_quantity(labels, quantity=None):
        try:
            cols = [colnames.index(label) for label in labels]
            if quantity:
                return convert(data[:, cols], quantity, units, "ASE")

            return data[:, cols]
        except ValueError:
            return None

    # slice data block into columns
    # + perform necessary conversions to ASE units
    positions = get_quantity(["x", "y", "z"], "distance")
    scaled_positions = get_quantity(["xs", "ys", "zs"])
    velocities = get_quantity(["vx", "vy", "vz"], "velocity")
    charges = get_quantity(["q"], "charge")
    forces = get_quantity(["fx", "fy", "fz"], "force")
    # !TODO: how need quaternions be converted?
    quaternions = get_quantity(["c_q[1]", "c_q[2]", "c_q[3]", "c_q[4]"])

    # convert cell
    cell = convert(cell, "distance", units, "ASE")
    celldisp = convert(celldisp, "distance", units, "ASE")
    if prismobj:
        celldisp = prismobj.vector_to_ase(celldisp)
        cell = prismobj.update_cell(cell)

    if quaternions:
        out_atoms = Quaternions(
            symbols=types,
            positions=positions,
            cell=cell,
            celldisp=celldisp,
            pbc=pbc,
            quaternions=quaternions,
        )
    elif positions is not None:
        # reverse coordinations transform to lammps system
        # (for all vectors = pos, vel, force)
        if prismobj:
            positions = prismobj.vector_to_ase(positions, wrap=True)

        out_atoms = atomsobj(
            symbols=types, positions=positions, pbc=pbc, celldisp=celldisp, cell=cell
        )
    elif scaled_positions is not None:
        out_atoms = atomsobj(
            symbols=types,
            scaled_positions=scaled_positions,
            pbc=pbc,
            celldisp=celldisp,
            cell=cell,
        )

    if velocities is not None:
        if prismobj:
            velocities = prismobj.vector_to_ase(velocities)
        out_atoms.set_velocities(velocities)
    if charges is not None:
        out_atoms.set_initial_charges(charges)
    if forces is not None:
        if prismobj:
            forces = prismobj.vector_to_ase(forces)
        # !TODO: use another calculator if available (or move forces
        #        to atoms.property) (other problem: synchronizing
        #        parallel runs)
        calculator = SinglePointCalculator(out_atoms, energy=0.0, forces=forces)
        out_atoms.set_calculator(calculator)

    return out_atoms


def construct_cell(diagdisp, offdiag):
    """Help function to create an ASE-cell with displacement vector from
    the lammps coordination system parameters.

    :param diagdisp: cell dimension convoluted with the displacement vector
    :param offdiag: off-diagonal cell elements
    :returns: cell and cell displacement vector
    :rtype: tuple
    """
    xlo, xhi, ylo, yhi, zlo, zhi = diagdisp
    xy, xz, yz = offdiag

    # create ase-cell from lammps-box
    xhilo = (xhi - xlo) - abs(xy) - abs(xz)
    yhilo = (yhi - ylo) - abs(yz)
    zhilo = zhi - zlo
    celldispx = xlo - min(0, xy) - min(0, xz)
    celldispy = ylo - min(0, yz)
    celldispz = zlo
    cell = np.array([[xhilo, 0, 0], [xy, yhilo, 0], [xz, yz, zhilo]])
    celldisp = np.array([celldispx, celldispy, celldispz])

    return cell, celldisp


def get_max_index(index):
    if np.isscalar(index):
        return index
    elif isinstance(index, slice):
        return index.stop if (index.stop is not None) else float("inf")


def read_lammps_dump_text(fileobj, index=-1, **kwargs):
    """Process cleartext lammps dumpfiles

    :param fileobj: filestream providing the trajectory data
    :param index: integer or slice object (default: get the last timestep)
    :returns: list of Atoms objects
    :rtype: list
    """
    # Load all dumped timesteps into memory simultaneously
    lines = deque(fileobj.readlines())

    index_end = get_max_index(index)

    n_atoms = 0
    images = []

    while len(lines) > n_atoms:
        line = lines.popleft()

        if "ITEM: TIMESTEP" in line:
            n_atoms = 0
            line = lines.popleft()
            # !TODO: pyflakes complains about this line -> do something
            # ntimestep = int(line.split()[0])  # NOQA

        if "ITEM: NUMBER OF ATOMS" in line:
            line = lines.popleft()
            n_atoms = int(line.split()[0])

        if "ITEM: BOX BOUNDS" in line:
            # save labels behind "ITEM: BOX BOUNDS" in triclinic case
            # (>=lammps-7Jul09)
            # !TODO: handle periodic boundary conditions in tilt_items
            tilt_items = line.split()[3:]
            celldatarows = [lines.popleft() for _ in range(3)]
            celldata = np.loadtxt(celldatarows)
            diagdisp = celldata[:, :2].reshape(6, 1).flatten()

            # determine cell tilt (triclinic case!)
            if len(celldata[0]) > 2:
                # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS"
                # to assign tilt (vector) elements ...
                offdiag = celldata[:, 2]
                # ... otherwise assume default order in 3rd column
                # (if the latter was present)
                if len(tilt_items) >= 3:
                    sort_index = [tilt_items.index(i) for i in ["xy", "xz", "yz"]]
                    offdiag = offdiag[sort_index]
            else:
                offdiag = (0.0,) * 3

            cell, celldisp = construct_cell(diagdisp, offdiag)

            # Handle pbc conditions
            if len(tilt_items) > 3:
                pbc = ["p" in d.lower() for d in tilt_items[3:]]
            else:
                pbc = (False,) * 3

        if "ITEM: ATOMS" in line:
            colnames = line.split()[2:]
            datarows = [lines.popleft() for _ in range(n_atoms)]
            data = np.loadtxt(datarows)
            out_atoms = lammps_data_to_ase_atoms(
                data=data,
                colnames=colnames,
                cell=cell,
                celldisp=celldisp,
                atomsobj=Atoms,
                pbc=pbc,
                **kwargs
            )
            images.append(out_atoms)

        if len(images) > index_end >= 0:
            break

    return images[index]


def read_lammps_dump_binary(
    fileobj, index=-1, colnames=None, intformat="SMALLBIG", **kwargs
):
    """Read binary dump-files (after binary2txt.cpp from lammps/tools)

    :param fileobj: file-stream containing the binary lammps data
    :param index: integer or slice object (default: get the last timestep)
    :param colnames: data is columns and identified by a header
    :param intformat: lammps support different integer size.  Parameter set \
    at compile-time and can unfortunately not derived from data file
    :returns: list of Atoms-objects
    :rtype: list
    """
    # depending on the chosen compilation flag lammps uses either normal
    # integers or long long for its id or timestep numbering
    # !TODO: tags are cast to double -> missing/double ids (add check?)
    tagformat, bigformat = dict(
        SMALLSMALL=("i", "i"), SMALLBIG=("i", "q"), BIGBIG=("q", "q")
    )[intformat]

    index_end = get_max_index(index)

    # Standard columns layout from lammpsrun
    if not colnames:
        colnames = ["id", "type", "x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz"]

    images = []

    # wrap struct.unpack to raise EOFError
    def read_variables(string):
        obj_len = struct.calcsize(string)
        data_obj = fileobj.read(obj_len)
        if obj_len != len(data_obj):
            raise EOFError
        return struct.unpack(string, data_obj)

    while True:
        try:
            # read header
            ntimestep, = read_variables("=" + bigformat)
            n_atoms, triclinic = read_variables("=" + bigformat + "i")
            boundary = read_variables("=6i")
            diagdisp = read_variables("=6d")
            if triclinic != 0:
                offdiag = read_variables("=3d")
            else:
                offdiag = (0.0,) * 3
            size_one, nchunk = read_variables("=2i")
            if len(colnames) != size_one:
                raise ValueError("Provided columns do not match binary file")

            # lammps cells/boxes can have different boundary conditions on each
            # sides (makes mainly sense for different non-periodic conditions
            # (e.g. [f]ixed and [s]hrink for a irradiation simulation))
            # periodic case: b 0 = 'p'
            # non-peridic cases 1: 'f', 2 : 's', 3: 'm'
            pbc = np.sum(np.array(boundary).reshape((3, 2)), axis=1) == 0

            cell, celldisp = construct_cell(diagdisp, offdiag)

            data = []
            for _ in range(nchunk):
                # number-of-data-entries
                n_data, = read_variables("=i")
                # retrieve per atom data
                data += read_variables("=" + str(n_data) + "d")
            data = np.array(data).reshape((-1, size_one))

            # map data-chunk to ase atoms
            out_atoms = lammps_data_to_ase_atoms(
                data=data,
                colnames=colnames,
                cell=cell,
                celldisp=celldisp,
                pbc=pbc,
                **kwargs
            )

            images.append(out_atoms)

            # stop if requested index has been found
            if len(images) > index_end >= 0:
                break

        except EOFError:
            break

    return images[index]
