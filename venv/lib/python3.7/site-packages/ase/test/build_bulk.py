from ase.data import chemical_symbols, reference_states
from ase.build import bulk


lat_map = dict(fcc='FCC',
               bcc='BCC',
               hcp='HEX',
               bct='BCT',
               diamond='FCC',
               #sc='CUB',
               #orthorhombic='ORC',
               rhombohedral='RHL')
lat_counts = {}


for Z, ref in enumerate(reference_states):
    if ref is None:
        continue

    structure = ref['symmetry']
    if structure not in lat_map:
        continue

    sym = chemical_symbols[Z]
    if sym in {'B', 'Se', 'Te'}:
        continue
    lat_counts.setdefault(structure, []).append(sym)

    atoms = bulk(sym)
    lat = atoms.cell.get_bravais_lattice()
    print(Z, atoms.symbols[0], structure, lat, atoms.cell.lengths())
    par1 = lat.tocell().niggli_reduce()[0].cellpar()
    par2 = atoms.cell.niggli_reduce()[0].cellpar()
    assert abs(par2 - par1).max() < 1e-10
    assert lat_map[structure] == lat.name

    if lat.name in ['RHL', 'BCT']:
        continue

    orth_atoms = bulk(sym, orthorhombic=True)
    orc_lat = orth_atoms.cell.get_bravais_lattice()
    angles = orc_lat.cellpar()[3:]
    assert abs(angles - 90).max() < 1e-10

    if lat.name in ['HEX', 'TET', 'ORC']:
        continue

    cub_atoms = bulk(sym, cubic=True)
    cub_lat = cub_atoms.cell.get_bravais_lattice()
    assert cub_lat.name == 'CUB', cub_lat


for key, val in lat_counts.items():
    print(key, len(val), ''.join(val))
