
import os
from ase.calculators.siesta.siesta import Siesta
from ase.constraints import FixAtoms, FixedLine, FixedPlane
from ase import Atoms


pseudo_path = 'pseudos'
if not os.path.exists(pseudo_path): os.makedirs(pseudo_path)

# Make dummy pseudopotentials.
for symbol in 'HCO':
    with open('{0}/{1}.lda.psf'.format(pseudo_path, symbol), 'w') as fd:
        fd.close()

atoms = Atoms('CO2', [(0.0, 0.0, 0.0), (-1.178, 0.0, 0.0), (1.178, 0.0, 0.0)])

c1 = FixAtoms(indices=[0])
c2 = FixedLine(1, [0.0, 1.0, 0.0])
c3 = FixedPlane(2, [1.0, 0.0, 0.0])

atoms.set_constraint([c1,c2,c3])

custom_dir = './dir1/'

# Test simple fdf-argument case.
siesta = Siesta(
    label=custom_dir + 'test_label',
    symlink_pseudos=False,
    atomic_coord_format='zmatrix',
    fdf_arguments={
        'MD.TypeOfRun': 'CG',
        'MD.NumCGsteps': 1000
        })

atoms.set_calculator(siesta)
siesta.write_input(atoms, properties=['energy'])

assert os.path.isfile(os.path.join(custom_dir, 'C.lda.1.psf'))
assert os.path.isfile(os.path.join(custom_dir, 'O.lda.2.psf'))

with open(os.path.join(custom_dir, 'test_label.fdf'), 'r') as fd: lines = fd.readlines()
lsl = [line.split() for line in lines]
assert ['cartesian'] in lsl
assert ['%block', 'Zmatrix'] in lsl
assert ['%endblock', 'Zmatrix'] in lsl
assert ['MD.TypeOfRun', 'CG'] in lsl
assert any([line.split()[4:9] == ['0', '0', '0', '1', 'C'] for line in lines])
assert any([line.split()[4:9] == ['0', '1', '0', '2', 'O'] for line in lines])
assert any([line.split()[4:9] == ['0', '1', '1', '3', 'O'] for line in lines])
