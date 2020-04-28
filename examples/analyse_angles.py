#!/usr/bin/env python3

'''
TODO: Description Needed

'''

def test_analyse_angles():

    from carmm.analyse.angles import analyse_all_angles

    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)
    #########

    analyse_all_angles(slab)
    # TODO: Add assertion test

test_analyse_angles()