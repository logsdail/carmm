from ase.calculators.test import FreeElectrons
from ase.lattice import all_variants
from ase.dft.band_structure import calculate_band_structure
from ase.utils import workdir
from ase import Atoms
import matplotlib.pyplot as plt

def test():
    ax = plt.gca()

    for i, lat in enumerate(all_variants()):
        if lat.ndim == 2:
            break
        xid = '{:02d}.{}'.format(i, lat.variant)
        path = lat.bandpath(density=10)
        path.write('path.{}.json'.format(xid))
        atoms = Atoms(cell=lat.tocell(), pbc=True)
        atoms.calc = FreeElectrons(nvalence=0, kpts=path.kpts)
        bs = calculate_band_structure(atoms, path)
        bs.write('bs.{}.json'.format(xid))
        bs.plot(ax=ax, emin=0, emax=20, filename='fig.{}.png'.format(xid))
        ax.clear()

with workdir('files', mkdir=True):
    test()
