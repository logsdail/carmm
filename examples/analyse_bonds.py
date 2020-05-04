#!/usr/bin/env python3

'''
This gives example usage when looking at how to get bond information from an Atoms object

This comes in useful when analysing input/output structures
'''

def test_analyse_bonds():

    from carmm.analyse.bonds import analyse_all_bonds

    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)
    #########

    analyse_all_bonds(slab)
    # TODO: Add assertion test

test_analyse_bonds()