'''
This example shows how to use the work flow functionality with CARMM to optimise
reactants/products and find transition states in an automated manner.
@Igor: Can you update with any further information?
'''

def test_run_workflows_ReactAims():
    from carmm.run.workflows.react import ReactAims
    from ase.build import molecule
    from carmm.analyse.forces import is_converged
    import os

    # Work in a dedicated folder
    parent_dir = os.getcwd()
    os.makedirs("data/react", exist_ok=True)
    os.chdir("data/react")

    # Determine calculation input settings
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

    # Create a reaction pathway
    atoms[1].x += 8
    transition_state = reactor.search_ts(atoms, model_optimised, 0.05, 0.03, input_check=0.01)
    activation_energy = transition_state.get_potential_energy()-model_optimised.get_potential_energy()

    assert 6.71774 == round(activation_energy, 5)

    # Calculate charges
    molecule_with_charges = reactor.get_mulliken_charges(model_optimised)
    assert molecule_with_charges[0].charge == 0.0

    '''The below uses the "dry_run" flag and uses a LennardJones calculator instead of FHI-aims to test code in CI'''
    reactor = ReactAims(params, basis_set, hpc, dry_run=True, filename="O2")

    '''Prepare TS search input'''
    from ase.build import molecule
    initial = molecule("O2")
    final = molecule("O2")
    final[1].x += 8

    # Provide one optimised image and one not converged - check using input_check
    # Then search for Transition State
    initial = reactor.aims_optimise(initial, 0.01)[0]

    '''Parallel task-farmed FHI-aims setup. The total number of images is n + 2 (middle images + input)
    Make sure total number of nodes requested in job submission is equal to nodes_per_instance * n. 
    E.g. for a band of 9 images and 1 node used per FHI-aims instance request 7 nodes in job submission'''

    reactor.nodes_per_instance = 1
    transition_state = reactor.search_ts_taskfarm(initial, final, 0.05, n=7, input_check=0.01)

    # Optimise a bulk geometry using stress tensor calculations and ExpCellFilter
    from ase.build import bulk
    reactor.nodes_per_instance = None
    reactor.filename = "Au"
    Au_bulk = bulk("Au")
    reactor.aims_optimise(Au_bulk, 0.01, relax_unit_cell=True)

    # Return to parent directory
    os.chdir(parent_dir)

test_run_workflows_ReactAims()
