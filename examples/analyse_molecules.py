def test_molecules():
    '''
    Test function for number of molecules in model

    '''
    from ase.visualize import view
    from carmm.analyse.molecules import calculate_molecules
    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    # Calculate neighbours
    molecules = calculate_molecules(slab)

    # View molecules
    A_mol = atoms[molecules[0]]
    B_mol = atoms[molecules[1]]
    #view(A_mol)
    #view(B_mol)
    #view(slab)
    # Verify results
    assert(len(molecules) == 2)




test_molecules()