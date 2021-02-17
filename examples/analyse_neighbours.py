def test_neighbours():
    '''Test neighgbour function

    TODO: function is a bit broken, running the following gives output 'None' '''
    from carmm.analyse.neighbours import neighbours
    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)
    neighbour_atoms = neighbours(slab,0,1)
    print(neighbour_atoms)

test_neighbours()