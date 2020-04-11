from ase.build import bulk
from ase.calculators.siesta import Siesta
from ase.io.siesta import _read_fdf_lines

atoms = bulk('Ti')
calc = Siesta()
atoms.calc = calc
calc.write_input(atoms, properties=['energy'])
# Should produce siesta.fdf but really we should be more explicit

fname = 'siesta.fdf'

with open(fname) as fd:
    thing = _read_fdf_lines(fd)
print(thing)

assert thing[0].split() == ['SystemName', 'siesta']
