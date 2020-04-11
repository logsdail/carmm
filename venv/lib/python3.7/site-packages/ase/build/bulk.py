from math import sqrt

from ase.atoms import Atoms
from ase.symbols import string2symbols
from ase.data import reference_states, atomic_numbers, chemical_symbols
from ase.utils import plural


def incompatible_cell(*, want, have):
    return RuntimeError('Cannot create {} cell for {} structure'
                        .format(want, have))


def bulk(name, crystalstructure=None, a=None, b=None, c=None, *, alpha=None,
         covera=None, u=None, orthorhombic=False, cubic=False,
         basis=None):
    """Creating bulk systems.

    Crystal structure and lattice constant(s) will be guessed if not
    provided.

    name: str
        Chemical symbol or symbols as in 'MgO' or 'NaCl'.
    crystalstructure: str
        Must be one of sc, fcc, bcc, hcp, diamond, zincblende,
        rocksalt, cesiumchloride, fluorite or wurtzite.
    a: float
        Lattice constant.
    b: float
        Lattice constant.  If only a and b is given, b will be interpreted
        as c instead.
    c: float
        Lattice constant.
    alpha: float
        Angle in degrees for rhombohedral lattice.
    covera: float
        c/a ratio used for hcp.  Default is ideal ratio: sqrt(8/3).
    u: float
        Internal coordinate for Wurtzite structure.
    orthorhombic: bool
        Construct orthorhombic unit cell instead of primitive cell
        which is the default.
    cubic: bool
        Construct cubic unit cell if possible.
    """

    if c is None and b is not None:
        # If user passes (a, b) positionally, we want it as (a, c) instead:
        c, b = b, c

    if covera is not None and c is not None:
        raise ValueError("Don't specify both c and c/a!")

    xref = None
    ref = {}

    if name in chemical_symbols:
        Z = atomic_numbers[name]
        ref = reference_states[Z]
        if ref is not None:
            xref = ref['symmetry']

            # If user did not specify crystal structure, and no basis
            # is given, and the reference state says we need one, but
            # does not have one, then we can't proceed.
            if (crystalstructure is None and basis is None
                and 'basis' in ref and ref['basis'] is None):
                # XXX This is getting much too complicated, we need to split
                # this function up.  A lot.
                raise RuntimeError('This structure requires an atomic basis')

        if ref is None:
            ref = {}  # easier to 'get' things from empty dictionary than None

        if xref == 'cubic':
            # P and Mn are listed as 'cubic' but the lattice constants
            # are 7 and 9.  They must be something other than simple cubic
            # then. We used to just return the cubic one but that must
            # have been wrong somehow.  --askhl
            raise RuntimeError('Only simple cubic ("sc") supported')

    # Mapping of name to number of atoms in primitive cell.
    structures = {'sc': 1, 'fcc': 1, 'bcc': 1,
                  'tetragonal': 1,
                  'bct': 1,
                  'hcp': 1,
                  'rhombohedral': 1,
                  'orthorhombic': 1,
                  'mcl': 1,
                  'diamond': 1,
                  'zincblende': 2, 'rocksalt': 2, 'cesiumchloride': 2,
                  'fluorite': 3, 'wurtzite': 2}

    if crystalstructure is None:
        crystalstructure = xref
        if crystalstructure not in structures:
            raise ValueError('No suitable reference data for bulk {}.'
                             '  Reference data: {}'
                             .format(name, ref))

    if crystalstructure not in structures:
        raise ValueError('Unknown structure: {}.'
                         .format(crystalstructure))

    # Check name:
    natoms = len(string2symbols(name))
    natoms0 = structures[crystalstructure]
    if natoms != natoms0:
        raise ValueError('Please specify {} for {} and not {}'
                         .format(plural(natoms0, 'atom'),
                                 crystalstructure, natoms))

    if alpha is None:
        alpha = ref.get('alpha')

    if a is None:
        if xref != crystalstructure:
            raise ValueError('You need to specify the lattice constant.')
        try:
            a = ref['a']
        except KeyError:
            raise KeyError('No reference lattice parameter "a" for "{}"'
                           .format(name))

    if b is None:
        bovera = ref.get('b/a')
        if bovera is not None and a is not None:
            b = bovera * a

    if crystalstructure in ['hcp', 'wurtzite']:
        if cubic:
            raise incompatible_cell(want='cubic', have=crystalstructure)

        if c is not None:
            covera = c / a
        elif covera is None:
            if xref == crystalstructure:
                covera = ref['c/a']
            else:
                covera = sqrt(8 / 3)

    if covera is None:
        covera = ref.get('c/a')
        if c is None and covera is not None:
            c = covera * a

    if orthorhombic and crystalstructure not in ['sc', 'tetragonal',
                                                 'orthorhombic']:
        return _orthorhombic_bulk(name, crystalstructure, a, covera, u)

    if cubic and crystalstructure in ['bcc', 'cesiumchloride']:
        return _orthorhombic_bulk(name, crystalstructure, a, covera)

    if cubic and crystalstructure != 'sc':
        return _cubic_bulk(name, crystalstructure, a)

    if crystalstructure == 'sc':
        atoms = Atoms(name, cell=(a, a, a), pbc=True)
    elif crystalstructure == 'fcc':
        b = a / 2
        atoms = Atoms(name, cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=True)
    elif crystalstructure == 'bcc':
        b = a / 2
        atoms = Atoms(name, cell=[(-b, b, b), (b, -b, b), (b, b, -b)],
                      pbc=True)
    elif crystalstructure == 'hcp':
        atoms = Atoms(2 * name,
                      scaled_positions=[(0, 0, 0),
                                        (1 / 3, 2 / 3, 0.5)],
                      cell=[(a, 0, 0),
                            (-a / 2, a * sqrt(3) / 2, 0),
                            (0, 0, covera * a)],
                      pbc=True)
    elif crystalstructure == 'diamond':
        atoms = bulk(2 * name, 'zincblende', a)
    elif crystalstructure == 'zincblende':
        s1, s2 = string2symbols(name)
        atoms = bulk(s1, 'fcc', a) + bulk(s2, 'fcc', a)
        atoms.positions[1] += a / 4
    elif crystalstructure == 'rocksalt':
        s1, s2 = string2symbols(name)
        atoms = bulk(s1, 'fcc', a) + bulk(s2, 'fcc', a)
        atoms.positions[1, 0] += a / 2
    elif crystalstructure == 'cesiumchloride':
        s1, s2 = string2symbols(name)
        atoms = bulk(s1, 'sc', a) + bulk(s2, 'sc', a)
        atoms.positions[1, :] += a / 2
    elif crystalstructure == 'fluorite':
        s1, s2, s3 = string2symbols(name)
        atoms = bulk(s1, 'fcc', a) + bulk(s2, 'fcc', a) + bulk(s3, 'fcc', a)
        atoms.positions[1, :] += a / 4
        atoms.positions[2, :] += a * 3 / 4
    elif crystalstructure == 'wurtzite':
        u = u or 0.25 + 1 / 3 / covera**2
        atoms = Atoms(2 * name,
                      scaled_positions=[(0, 0, 0),
                                        (1 / 3, 2 / 3, 0.5 - u),
                                        (1 / 3, 2 / 3, 0.5),
                                        (0, 0, 1 - u)],
                      cell=[(a, 0, 0),
                            (-a / 2, a * sqrt(3) / 2, 0),
                            (0, 0, a * covera)],
                      pbc=True)
    elif crystalstructure == 'bct':
        from ase.lattice import BCT
        if basis is None:
            basis = ref.get('basis')
        if basis is not None:
            natoms = len(basis)
        lat = BCT(a=a, c=c)
        atoms = Atoms([name] * natoms, cell=lat.tocell(), pbc=True,
                      scaled_positions=basis)
    elif crystalstructure == 'rhombohedral':
        atoms = _build_rhl(name, a, alpha, basis)
    elif crystalstructure == 'orthorhombic':
        atoms = Atoms(name, cell=[a, b, c], pbc=True)
    else:
        raise ValueError('Unknown crystal structure: ' + crystalstructure)

    if orthorhombic:
        assert atoms.cell.orthorhombic
    if cubic:
        assert abs(atoms.cell.angles() - 90).all() < 1e-10
    return atoms


def _build_rhl(name, a, alpha, basis):
    from ase.lattice import RHL
    lat = RHL(a, alpha)
    cell = lat.tocell()
    if basis is None:
        # RHL: Given by A&M as scaled coordinates "x" of cell.sum(0):
        basis_x = reference_states[atomic_numbers[name]]['basis_x']
        basis = basis_x[:, None].repeat(3, axis=1)
    natoms = len(basis)
    return Atoms([name] * natoms, cell=cell, scaled_positions=basis, pbc=True)


def _orthorhombic_bulk(name, crystalstructure, a, covera=None, u=None):
    if crystalstructure == 'fcc':
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    elif crystalstructure == 'bcc':
        atoms = Atoms(2 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    elif crystalstructure == 'hcp':
        atoms = Atoms(4 * name,
                      cell=(a, a * sqrt(3), covera * a),
                      scaled_positions=[(0, 0, 0),
                                        (0.5, 0.5, 0),
                                        (0.5, 1 / 6, 0.5),
                                        (0, 2 / 3, 0.5)],
                      pbc=True)
    elif crystalstructure == 'diamond':
        atoms = _orthorhombic_bulk(2 * name, 'zincblende', a)
    elif crystalstructure == 'zincblende':
        s1, s2 = string2symbols(name)
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0, 0.25),
                                        (0.5, 0.5, 0.5), (0, 0.5, 0.75)])
    elif crystalstructure == 'rocksalt':
        s1, s2 = string2symbols(name)
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0),
                                        (0.5, 0.5, 0.5), (0, 0, 0.5)])
    elif crystalstructure == 'cesiumchloride':
        atoms = Atoms(name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    elif crystalstructure == 'wurtzite':
        u = u or 0.25 + 1 / 3 / covera**2
        atoms = Atoms(4 * name,
                      cell=(a, a * 3**0.5, covera * a),
                      scaled_positions=[(0, 0, 0),
                                        (0, 1 / 3, 0.5 - u),
                                        (0, 1 / 3, 0.5),
                                        (0, 0, 1 - u),
                                        (0.5, 0.5, 0),
                                        (0.5, 5 / 6, 0.5 - u),
                                        (0.5, 5 / 6, 0.5),
                                        (0.5, 0.5, 1 - u)],
                      pbc=True)
    else:
        raise incompatible_cell(want='orthorhombic', have=crystalstructure)

    return atoms


def _cubic_bulk(name, crystalstructure, a):
    if crystalstructure == 'fcc':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0, 0.5, 0.5),
                                        (0.5, 0, 0.5), (0.5, 0.5, 0)])
    elif crystalstructure == 'diamond':
        atoms = _cubic_bulk(2 * name, 'zincblende', a)
    elif crystalstructure == 'zincblende':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.25, 0.25, 0.25),
                                        (0, 0.5, 0.5), (0.25, 0.75, 0.75),
                                        (0.5, 0, 0.5), (0.75, 0.25, 0.75),
                                        (0.5, 0.5, 0), (0.75, 0.75, 0.25)])
    elif crystalstructure == 'rocksalt':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0, 0),
                                        (0, 0.5, 0.5), (0.5, 0.5, 0.5),
                                        (0.5, 0, 0.5), (0, 0, 0.5),
                                        (0.5, 0.5, 0), (0, 0.5, 0)])
    else:
        raise incompatible_cell(want='cubic', have=crystalstructure)

    return atoms
