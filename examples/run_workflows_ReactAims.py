def test_react():
    from carmm.run.workflows.react import ReactAims
    from ase.build import molecule
    from carmm.analyse.forces import is_converged

    # Determine calculation input settings
    params = {"xc":"pbe"}
    basis_set = "light"
    hpc = "hawk"

    '''Prepare the Atoms object geometry'''
    atoms = molecule("H2O")

    '''Build the React_Aims object using the set of settings'''
    reactor = ReactAims(params, basis_set, hpc)

    '''Call relevant calculations'''
    model_optimised, model_postprocessed = reactor.aims_optimise(atoms, fmax=0.01, restart=True)

    # TODO: The ase.vibrations module was changed between 3.21.0 and 3.22.1, the below does not work in old version
    zero_point_energy = reactor.vibrate(atoms, indices =[atom.index for atom in atoms])

    # TODO: figure out a test for ts_search and get_mulliken_charge

    assert is_converged(reactor.model_optimised, 0.01), \
        "The structure saved in React_Aims is not converged"
    assert round(zero_point_energy, 3) == 0.574
    assert atoms.calc == None

test_react()