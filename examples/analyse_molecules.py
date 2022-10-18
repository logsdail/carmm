def test_molecules():
    '''
    Test function for number of molecules in model

    '''
    from ase.visualize import view
    from carmm.analyse.molecules import calculate_molecules, calculate_formula
    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    # Calculate neighbours
    molecules = calculate_molecules(slab)
    symbols = calculate_formula(slab)
    print(symbols)
    # View seperate molecules and recombine
    A_mol = slab[molecules[0]]
    B_mol = slab[molecules[1]]
    #view(A_mol)
    #view(B_mol)
    #view(slab)
    # Verify results
    assert(len(molecules) == 2)
    assert(symbols[0] == 'Au18' and symbols[1] == 'C1O2')




test_molecules()