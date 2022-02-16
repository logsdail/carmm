from carmm.run.workflows.react import React_Aims
from ase.io import read

# TODO: turn into an example for the regression test
'''Setup FHI-aims calculation parameters'''
params = {"xc":"pbe"}
dimensions = 0
basis_set = "light"
hpc = "hawk"


'''Prepare the Atoms object geometry'''
atoms = read()
#initial = read()
#final = read()
'''Run the ASE/FHI-aims BFGS force optimisation on the hpc facility'''
reactor_opt = React_Aims(params, basis_set, hpc, dimensions)
reactor_opt.aims_optimise(atoms,
                      fmax=0.01,
                      post_process=None,
                      relax_unit_cell=False,
                      internal= False,
                      restart= False)

'''Run the ASE/FHI-aims/CatLearn transition state search on the hpc facility'''
#reactor_ts = React_Aims(params, basis_set, hpc, dimensions)
#reactor_ts.search_ts(initial, final, 0.05, 0.03)

