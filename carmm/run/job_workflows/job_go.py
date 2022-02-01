from ase.io import read
from geometry_optimisation import aims_optimise
from ase.constraints import FixAtoms
from carmm.analyse.neighbours import neighbours

hpc = "archer2"

model = read()

"""
# Sequential constraints - neighbour shells
adsorbate_indices = [atom.index for atom in model if atom.symbol in ["C", "H", "O"]]
bottom_layer = FixAtoms([atom.index for atom in model if atom.tag > 4])
first_shell = FixAtoms(
    [atom.index for atom in model if atom.index not in neighbours(model, adsorbate_indices, 1)[0]])

# Preoptimisation with just adsorbate + 1st neighbours
model = aims_optimise(model, hpc, [first_shell, bottom_layer], fmax=0.03, tight=False, preopt=True)[0]

second_shell = FixAtoms(
    [atom.index for atom in model if atom.index not in neighbours(model, adsorbate_indices, 2)[0]])

# Preoptimisation with adsorbate and 2 shells of neighbours
model = aims_optimise(model, hpc, [second_shell, bottom_layer], fmax=0.03, tight=False, preopt=True)[0]

# Optimisation with a fully relaxed surface and adsorbate (bulk layers still constrained)
model = aims_optimise(model, hpc, [bottom_layer], fmax=0.01, dimensions=2, tight=True, preopt=False)[0]
"""

# Optimisation of a bulk with unit cell relaxation
model = aims_optimise(model, hpc, fmax=0.01, dimensions=3, sf=True)[0]

print()
print("Success.")