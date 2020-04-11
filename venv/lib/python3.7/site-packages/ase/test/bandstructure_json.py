from ase.build import bulk
from ase.dft.band_structure import calculate_band_structure, BandStructure
from ase.io.jsonio import read_json
from ase.calculators.test import FreeElectrons

atoms = bulk('Au')
lat = atoms.cell.get_bravais_lattice()
path = lat.bandpath(npoints=100)

atoms.calc = FreeElectrons()

bs = calculate_band_structure(atoms, path)
bs.write('bs.json')
bs.path.write('path.json')

bs1 = read_json('bs.json')
bs2 = BandStructure.read('bs.json')
path1 = read_json('path.json')
assert type(bs1) == type(bs)
assert type(bs2) == type(bs)
assert type(path1) == type(bs.path)
