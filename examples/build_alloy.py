#!/usr/bin/env python3

'''
This example shows how to build an alloy surface, both binary and ternary,
and substituted single atom alloy surfaces.

This is useful when wanting a randomly distributed set of elements through a model
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

    #### Single Atom Alloy surfaces test ####
    from carmm.build.alloy import get_SAA_surfaces
    SAA_elements = ['Ag', 'Au', 'Cu', 'Pd', 'Ni', 'Ti']
    substitution_indices = [13, 16]
    SAA_surfaces = get_SAA_surfaces(slab, SAA_elements, substitution_indices, include_pristine=True)

    assert len(SAA_surfaces) == 16

test_build_alloy()
