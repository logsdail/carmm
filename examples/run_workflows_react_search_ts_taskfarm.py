'''
This example shows how to use the work flow functionality with CARMM to optimise
reactants/products and find transition states in an automated manner using task-farming
in an ASE/i-Pi/FHI-aims setup
'''


def test_run_workflows_react_ts_taskfarm():
    from carmm.run.workflows.react import ReactAims
    from ase.build import bulk
    import os
    from ase.build import surface, add_adsorbate
    from ase.constraints import FixAtoms

    '''Work in a dedicated folder'''
    parent_dir = os.getcwd()
    os.makedirs("data/react", exist_ok=True)
    os.chdir("data/react")

    '''Determine calculation input settings'''
    params = {"xc": "pbe"}
    basis_set = "light"
    hpc = "archer2"

    '''Optimise the bulk metal using stress tensor calculations and ExpCellFilter to later cut a surface model'''
    reactor = ReactAims(params, basis_set, hpc)

    '''Cut a 2x2-Al(001) surface with 3 layers and an
     Au atom adsorbed in a hollow site:'''
    initial = surface(bulk("Pd"), (1, 0, 0), layers=1, vacuum=20)
    initial = initial.repeat((2, 2, 1))
    add_adsorbate(initial, 'Al', 2, position=initial[0].position[:2])

    '''Fix second and third layers:'''
    mask = [atom.tag > 0 for atom in initial]
    initial.set_constraint(FixAtoms(mask=mask))

    '''Adjust reciprocal space s2mpling for surface models'''
    reactor.params["k_grid"] = (2, 2, 1)
    '''Filename based on chemical formula by default'''
    reactor.filename = None

    '''TS input - Provide initial optimised image and one not converged - check using input_check. 
    Calculate the Transition State'''
    final = initial.copy()
    final[-1].x += final.get_cell()[0, 0] / 2

    '''Below is the the task-farmed FHI-aims setup. The total number of images is n + 2 (middle images + input)
    Make sure total number of nodes requested in job submission is equal to nodes_per_instance * n. 
    E.g. for a band of 7 images and 1 node used per FHI-aims instance request 5 nodes in job submission'''

    reactor.nodes_per_instance = 1
    TS_CINEB = reactor.search_ts_taskfarm(initial=initial,
                                          final=final,
                                          fmax=0.05,
                                          n=3,
                                          method='string',
                                          input_check=0.05,
                                          interpolation="idpp",
                                          max_steps=100,
                                          verbose=True)

    '''Return to parent directory'''
    os.chdir(parent_dir)


test_run_workflows_react_ts_taskfarm()

