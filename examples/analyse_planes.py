def test_analyse_planes():

    from carmm.analyse.planes import get_interplane_distances

    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)
    #########

    distances_sorted = sorted(get_interplane_distances(slab))
    # As the adsorbate is 3 Angstrom above the surface, this should be the shortest distance
    assert(1e-5 > abs(3.0 - distances_sorted[0]))

test_analyse_planes()
