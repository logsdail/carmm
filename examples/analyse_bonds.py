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
    # Deform the z-coordinate for C so the Au-C bond length is too short, to test deformation
    slab.positions[-3][2] -= 1.5
    #########

    abnormal_count, abnormal_bond_names = analyse_all_bonds(slab, abnormal=True)

    assert(len(abnormal_count) == 1)

test_analyse_bonds()