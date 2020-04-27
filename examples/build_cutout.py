#!/usr/bin/env python3

'''
TODO: Description Needed

'''

def test_cutout():

    from software.build.cutout import cutout_sphere

    #### Traditional ASE functionality #####
    from software.examples.data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)
    #########

    cutout = cutout_sphere(slab, 13)
    assert(len(cutout) == 12)

test_cutout()
