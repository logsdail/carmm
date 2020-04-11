
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.utils.forcecurve import force_curve
from ase.md import VelocityVerlet
from ase.units import fs
from ase.io import read

atoms = bulk('Au', cubic=True) * (2, 1, 1)
atoms.calc = EMT()
atoms.rattle(stdev=0.05)

md = VelocityVerlet(atoms, timestep=12.0 * fs, trajectory='tmp.traj')
md.run(steps=10)
images = read('tmp.traj', ':')
force_curve(images)

# import pylab as plt
# plt.show()
