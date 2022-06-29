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

    # View seperate molecules and recombine
    A_mol = slab[molecules[0]]
    B_mol = slab[molecules[1]]
    recombine = A_mol + B_mol
    #view(A_mol)
    #view(B_mol)
    #view(slab)
    #view(recombine)
    # Verify results
    assert(len(molecules) == 2)




test_molecules()