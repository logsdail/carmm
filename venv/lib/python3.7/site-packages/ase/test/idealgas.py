from ase.md import VelocityVerlet
from ase.build import bulk
from ase.units import kB
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.idealgas import IdealGas
import numpy as np

atoms = bulk('Kr').repeat((10,10,10))
assert len(atoms) == 1000

atoms.center(vacuum=100)
atoms.set_calculator(IdealGas())

T = 1000

MaxwellBoltzmannDistribution(atoms, T * kB)
print("Temperature: {} K".format(atoms.get_temperature()))

md = VelocityVerlet(atoms, timestep=0.1)
for i in range(5):
    md.run(5)
    s = atoms.get_stress(include_ideal_gas=True)
    p = -s[:3].sum()/3
    v = atoms.get_volume()
    N = len(atoms)
    T = atoms.get_temperature()
    print("pV = {}  NkT = {}".format(p*v, N*kB*T))
    assert np.fabs(p*v - N*kB*T) < 1e-6
    
