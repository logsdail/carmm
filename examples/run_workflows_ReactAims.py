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
    model_optimised, model_postprocessed = reactor.aims_optimise(atoms, fmax=0.01, restart=True)

    # TODO: The ase.vibrations module was changed between 3.21.0 and 3.22.1, the below does not work in old version
    zero_point_energy = reactor.vibrate(atoms, indices =[atom.index for atom in atoms])

    assert is_converged(reactor.model_optimised, 0.01), \
        "The structure saved in React_Aims is not converged"
    assert round(zero_point_energy, 3) == 0.275

    # Create a reaction pathway
    atoms[1].x += 8
    transition_state = reactor.search_ts(atoms, model_optimised, 0.05, 0.03, input_check=0.01)
    activation_energy = transition_state.get_potential_energy()-model_optimised.get_potential_energy()

    assert 6.71774 == round(activation_energy, 5)

    # Calculate charges
    molecule_with_charges = reactor.get_mulliken_charges(model_optimised)
    assert molecule_with_charges[0].charge == 0.0

    # Return to parent directory
    os.chdir(parent_dir)

test_run_workflows_ReactAims()
