from carmm.run.workflows.react import React_Aims
from ase.build import molecule

# TODO: turn into an example for the regression test
'''Setup FHI-aims calculation parameters'''
params = {"xc":"pbe"}
dimensions = 0
basis_set = "light"
hpc = "hawk"
reactor = React_Aims(params, basis_set, hpc, dimensions)

'''Prepare the Atoms object geometry'''
atoms = molecule("H2O")

'''Run the ASE/FHI-aims BFGS force optimisation on the hpc facility'''
reactor.aims_optimise(atoms,
                      fmax=0.01,
                      post_process=None,
                      relax_unit_cell=False,
                      internal= False,
                      restart= False)

print(reactor.model_optimised.get_potential_energy())
print(reactor.model_post_processed.get_potential_energy())
