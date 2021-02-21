def test_neighbours():
    '''
    Test neighbour function

    TODO: Check the given neighbours are correct!
    '''
    from carmm.analyse.neighbours import neighbours
    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    # Calculate neighbours
    neighbour_atoms = neighbours(slab,0,1)

    # Verify results
    assert(neighbour_atoms == [0, 1, 3, 9, 12])

test_neighbours()