
from ase.calculators.octopus import Octopus
from ase.calculators.interfacechecker import check_interface
from ase.build import molecule

system = molecule('H2')
system.center(vacuum=2.0)

directory = 'ink'

calc0 = Octopus(directory=directory,
                check_keywords=False,
                FromScratch=True,
                stdout="'stdout.txt'",
                stderr="'stderr.txt'",
                Spacing='0.25 * Angstrom',
                OutputFormat='cube + xcrysden')

system.set_calculator(calc0)
system.get_potential_energy()

# Must make one test with well-defined cell and one without.

calc1 = Octopus(directory)
system = calc1.get_atoms()

E = system.get_potential_energy()
print('energy', E)

errs = check_interface(calc1)
# view(system)

atoms = Octopus.read_atoms(directory)
errs = check_interface(atoms.calc)

changes = calc1.check_state(atoms)
print('changes', changes)
assert len(changes) == 0
