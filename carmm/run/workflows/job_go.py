from ase.io import read
from carmm.run.workflows.geometry_optimisation import aims_optimise
from carmm.run.aims_calculator import get_k_grid

# TODO: refactor into an example that can go into the regression test

'''Provide the name of the hpc facility you are working on.'''
hpc = "hawk"

'''Provide your structure in the form of an Atoms object, e.g. by reading a file.'''
model = read("input.traj")

'''Periodicity - 0 for gas, 2 for surface, 3 for bulk.'''
dimensions = 2

'''Define your FHI-aims parameters'''
params = {
         "xc_pre":['pbe', '10'],
         "xc":'libxc MGGA_X_MBEEF+GGA_C_PBE_SOL',
         "spin":'none',
         "relativistic":('atomic_zora','scalar'),
         #"force_correction":'True',
         #"sc_accuracy_etot":1e-6,
         #"sc_accuracy_eev":1e-3,
         #"sc_accuracy_rho":1e-6,
         #"sc_accuracy_forces":1e-4,
        }

'''Define k_grid for any periodic calculation, sampling density is in k-points per Å^-1'''
k_sampling_density = 0.018
if dimensions > 0:
    params["k_grid"] = get_k_grid(model, k_sampling_density, dimensions=dimensions, verbose=True)

'''Geometry optimisation with a smaller basis set and a single point calculation with larger basis set'''
'''The basis sets need to be adjusted to user's preference inn geometry_optimisation'''
model, model_tight = aims_optimise(params,
                                   model,
                                   dimensions,
                                   hpc,
                                   fmax=0.01,
                                   tight=True,
                                   sf= False,
                                   internal=False,
                                   restart=False
                                   )

print()
print("Success.")

'''Optimisation of a bulk with unit cell relaxation'''
# model = aims_optimise(params, model, hpc, fmax=0.01, dimensions=3, sf=True, restart=False)[0]

"""
'''Sequential optimisation setup'''
# TODO: This section must be separated
from ase.constraints import FixAtoms
from carmm.analyse.neighbours import neighbours

# Preoptimisation with just adsorbate + 1st neighbours
model = aims_optimise(model, hpc, [first_shell, bottom_layer], fmax=0.03, tight=False, preopt=True)[0]

second_shell = FixAtoms(
    [atom.index for atom in model if atom.index not in neighbours(model, adsorbate_indices, 2)[0]])

# Preoptimisation with adsorbate and 2 shells of neighbours
model = aims_optimise(model, hpc, [second_shell, bottom_layer], fmax=0.03, tight=False, preopt=True)[0]

# Optimisation with a fully relaxed surface and adsorbate (bulk layers still constrained)
model = aims_optimise(model, hpc, [bottom_layer], fmax=0.01, dimensions=2, tight=True, preopt=False)[0]
"""