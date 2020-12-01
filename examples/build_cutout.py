#!/usr/bin/env python3

'''
This example shows to to build a model where you "cutout" a spherical model from a slab with adsorbate

This is really useful when wanting to do BSSE calculations (as you can't do them with periodic models)
'''

def test_build_cutout():

    from carmm.build.cutout import cutout_sphere

    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)
    #########

    cutout = cutout_sphere(slab, 13)
    assert(len(cutout) == 12)

test_build_cutout()
