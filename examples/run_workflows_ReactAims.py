'''
This example shows how to use the work flow functionality with CARMM to optimise
reactants/products and find transition states in an automated manner.
@Igor: Can you update with any further information?
'''

def test_run_workflows_ReactAims():
    from carmm.run.workflows.react import ReactAims
    from ase.build import molecule, bulk
    from carmm.analyse.forces import is_converged
    import os
    from ase.build import surface, add_adsorbate
    from ase.constraints import FixAtoms

    '''Work in a dedicated folder'''
    parent_dir = os.getcwd()
    os.makedirs("data/react", exist_ok=True)
    os.chdir("data/react")

    '''Determine calculation input settings'''
    params = {"xc":"pbe"}
    basis_set = "light"
    hpc = "hawk"

    '''Prepare the Atoms object geometry'''
    atoms = molecule("H2")

    '''Build the React_Aims object using the set of settings'''
    reactor = ReactAims(params, basis_set, hpc)

    '''Call relevant calculations'''
    '''The below has been previously calculated and data is retrieved from saved trajectories'''
    model_optimised, model_postprocessed = reactor.aims_optimise(atoms, fmax=0.01, restart=True)
    zero_point_energy = reactor.vibrate(atoms, indices =[atom.index for atom in atoms])

    assert is_converged(reactor.model_optimised, 0.01), \
    '''The structure saved in React_Aims is not converged'''
    assert round(zero_point_energy, 3) == 0.275

    '''Create a reaction pathway'''
    atoms[1].x += 8
    transition_state = reactor.search_ts(atoms, model_optimised, 0.05, 0.03, input_check=0.01)
    activation_energy = transition_state.get_potential_energy()-model_optimised.get_potential_energy()

    assert 6.71774 == round(activation_energy, 5)

    '''Calculate charges'''
    molecule_with_charges = reactor.get_mulliken_charges(model_optimised)
    assert molecule_with_charges[0].charge == 0.0

    '''The below uses the "dry_run" flag and uses an EMT calculator instead of FHI-aims to test code in CI'''
    # TODO: Add relevant assertion statements below
    reactor = ReactAims(params, basis_set, hpc, dry_run=True, filename="O2")

    '''Optimise the bulk metal using stress tensor calculations and ExpCellFilter to later cut a surface model'''
    reactor.filename = "Al"
    reactor.params["k_grid"] = (8, 8, 8)
    Al_bulk = bulk("Al")

    '''Calculate the optimal unit cell and post process the calculation with a larger "tight" basis set'''
    light, tight = reactor.aims_optimise(Al_bulk, 0.01, relax_unit_cell=True, post_process="tight")

    '''Cut a 2x2-Al(001) surface with 3 layers and an
     Au atom adsorbed in a hollow site:'''
    initial = surface(Al_bulk, (1, 0, 0), layers=1, vacuum=10)
    initial = initial.repeat((2,1,1))
    add_adsorbate(initial, 'H', 1.7, position=initial[0].position[:2])


    '''Fix second and third layers:'''
    mask = [atom.tag > 0 for atom in initial]
    initial.set_constraint(FixAtoms(mask=mask))

    '''Adjust reciprocal space sampling for surface models'''
    reactor.params["k_grid"] = (4, 4, 1)
    '''Filename based on chemical formula by default'''
    reactor.filename = None

    '''TS input - Provide initial optimised image and one not converged - check using input_check. 
    Calculate the Transition State using ML-NEB'''
    initial = reactor.aims_optimise(initial, 0.01)[0]
    final = initial.copy()
    final[-1].x += final.get_cell()[0, 0] / 2

    TS_MLNEB = reactor.search_ts(initial, final, 0.05, 0.03, n=7, input_check=0.01, restart=True)

    '''Calculate the Transition State using AIDNEB from ase-gpatom package'''
    '''WARNING DO NOT USE WITH FHI-AIMS  - requires modified gpatom source code, issue opened on ase-gpatom GitHub'''

    TS_AIDNEB = reactor.search_ts_aidneb(initial, final, 0.05, 0.03, n=7, input_check=0.01, restart=True)

    '''Below is the the task-farmed FHI-aims setup. The total number of images is n + 2 (middle images + input)
    Make sure total number of nodes requested in job submission is equal to nodes_per_instance * n. 
    E.g. for a band of 7 images and 1 node used per FHI-aims instance request 5 nodes in job submission'''

    reactor.nodes_per_instance = 1
    TS_CINEB = reactor.search_ts_taskfarm(initial, final, 0.05, n=5, input_check=0.01, max_steps=30)

    '''Return to parent directory'''
    os.chdir(parent_dir)


test_run_workflows_ReactAims()

