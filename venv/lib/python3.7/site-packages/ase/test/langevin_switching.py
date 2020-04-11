import numpy as np
from ase.build import bulk
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.harmonic import SpringCalculator
from ase.md.switch_langevin import SwitchLangevin


# params
size = 6
T = 300
n_steps = 500
k1 = 2.0
k2 = 4.0
dt = 10

# for reproducibility
np.random.seed(42)

# setup atoms and calculators
atoms = bulk('Al').repeat(size)
calc1 = SpringCalculator(atoms.positions, k1)
calc2 = SpringCalculator(atoms.positions, k2)

# theoretical diff
n_atoms = len(atoms)
calc1.atoms = atoms
calc2.atoms = atoms
F1 = calc1.get_free_energy(T) / n_atoms
F2 = calc2.get_free_energy(T) / n_atoms
dF_theory = F2 - F1

# switch_forward
dyn_forward = SwitchLangevin(atoms, calc1, calc2, dt * units.fs, T * units.kB, 0.01, n_steps, n_steps)
MaxwellBoltzmannDistribution(atoms, 2 * T * units.kB)
dyn_forward.run()
dF_forward = dyn_forward.get_free_energy_difference() / len(atoms)

# switch_backwards
dyn_backward = SwitchLangevin(atoms, calc2, calc1, dt * units.fs, T * units.kB, 0.01, n_steps, n_steps)
MaxwellBoltzmannDistribution(atoms, 2 * T * units.kB)
dyn_backward.run()
dF_backward = -dyn_backward.get_free_energy_difference() / len(atoms)

# summary
dF_switch = (dF_forward + dF_backward) / 2.0
error = dF_switch - dF_theory

# print('delta_F analytical: {:12.6f} eV/atom'.format(dF_theory))
# print('delta_F forward:    {:12.6f} eV/atom'.format(dF_forward))
# print('delta_F backward:   {:12.6f} eV/atom'.format(dF_backward))
# print('delta_F average:    {:12.6f} eV/atom'.format(dF_switch))
# print('delta_F error:      {:12.6f} eV/atom'.format(error))
assert abs(error) < 1e-3
