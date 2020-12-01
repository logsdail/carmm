#!/usr/bin/env python3

'''
This example shows how to get a range of different surfaces for a material

This is useful when making surfaces for adsorption chemistry

TODO: Link with pymatgen, as that has extensive functionality?
'''

def test_build_surface():

    #### Functionality to create all surfaces ####
    from carmm.build.facets import generate
    facets, slabs = generate('Au', save=True)

    #### Assertion test for save ####
    import os
    assert(os.path.exists('111.in'))
    #########

    #### Assertion tests ####
    assert(len(slabs) == 1)
    assert(slabs[0].get_volume() - 644.5655232 < 1e-6)

# Run the example
test_build_surface()