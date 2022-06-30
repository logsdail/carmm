#!/usr/bin/env python3

'''
This short example shows how to get all the angluar information from an atoms object

This comes in useful when analysing input/output structures
'''

def test_analyse_angles():

    from carmm.analyse.angles import analyse_all_angles

    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)
    #########

    elements, indices, angles = analyse_all_angles(slab, verbose=True)

    index_au = elements.index(('Au', 'Au', 'Au'))
    assert(len(elements) == 2)
    assert(len(indices) == 2 and len(indices[index_au]) == 648)
    assert(len(angles) == 2 and len(angles[index_au]) == 648)

test_analyse_angles()

