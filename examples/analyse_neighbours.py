def test_neighbours():
    '''
    Test neighbour function
    '''
    from carmm.analyse.neighbours import neighbours
    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    from ase.visualize import view
    view(slab)

    # Calculate neighbours
    all_neighbour_atoms, shell_list = neighbours(slab, [13], 1)

    # Verify results
    assert(all_neighbour_atoms == [1, 2, 4, 10, 11, 12, 13, 14, 15, 16])
    assert(shell_list == [[13], [1, 2, 4, 10, 11, 12, 14, 15, 16]])
test_neighbours()