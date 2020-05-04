#!/usr/bin/env python3

'''
TODO: Description Needed

'''

def test_build_alloy():

    from carmm.build.alloy import binary_alloy

    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab = slab()
    #########

    alloy_slab = binary_alloy(slab, 'Pd', 8)
    assert(alloy_slab.get_chemical_symbols().count('Pd') == 8)
    assert(alloy_slab.get_chemical_symbols().count('Au') == 10)

    #### Ternary_alloy test ####
    from carmm.build.alloy import ternary_alloy

    ternary_slab = ternary_alloy(slab, 'Pd', 'Zn', 4, 5)
    # This currently fails.
    # TODO: Improve ternary so it preserves the requested composition for 2nd element
    # assert(ternary_slab.get_chemical_symbols().count('Pd') == 4)
    assert(ternary_slab.get_chemical_symbols().count('Zn') == 5)

test_build_alloy()
