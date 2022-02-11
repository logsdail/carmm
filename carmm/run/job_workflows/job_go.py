from ase.io import read
from geometry_optimisation import aims_optimise, get_k_grid
from ase.constraints import FixAtoms
from carmm.analyse.neighbours import neighbours

hpc = "hawk"
model = read("input.traj")

# periodicity - 0 for gas, 2 for surface, 3 for bulk
dimensions = 2

# Define your FHI-aims parameters
params = {
         "xc_pre":['pbe', '10'],
         "xc":'libxc MGGA_X_MBEEF+GGA_C_PBE_SOL',
         "spin":'none',
         "relativistic":('atomic_zora','scalar'),
         "force_correction":'True',
         #"sc_accuracy_etot":1e-6,
         #"sc_accuracy_eev":1e-3,
         #"sc_accuracy_rho":1e-6,
         #"sc_accuracy_forces":1e-4,
         #"final_forces_cleaned":'true',
        }

# define k_grid for any periodic calculation
k_sampling_density = 0.018 # k-points per Å^-1
if dimensions > 0:
    params["k_grid"] = get_k_grid(model, k_sampling_density, dimensions=dimensions, verbose=True)

# Optimisation of a bulk with unit cell relaxation
#model = aims_optimise(params, model, hpc, fmax=0.01, dimensions=3, sf=True, restart=False)[0]

# Generic optimisation
model = aims_optimise(params, model, hpc, fmax=0.01, dimensions=dimensions)[0]

print()
print("Success.")




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