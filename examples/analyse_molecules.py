def test_molecules():
    '''
    Test function for number of molecules in model

    '''
    from carmm.analyse.molecules import calculate_molecules
    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    # Calculate neighbours
    molecules = calculate_molecules(slab)

    # Verify results
    assert(len(molecules) == 2)

test_molecules()