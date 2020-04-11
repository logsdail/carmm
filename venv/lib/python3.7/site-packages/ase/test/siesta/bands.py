from ase.build import bulk
from ase.calculators.siesta import Siesta
from ase.utils import workdir
from ase.dft.band_structure import calculate_band_structure

atoms = bulk('Si')
with workdir('files', mkdir=True):
    path = atoms.cell.bandpath('GXWK', density=10)
    atoms.calc = Siesta(kpts=[2, 2, 2])
    bs = calculate_band_structure(atoms, path)
    print(bs)
    bs.write('bs.json')
