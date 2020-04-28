#!/usr/bin/env python3

'''
TODO: Description Needed

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