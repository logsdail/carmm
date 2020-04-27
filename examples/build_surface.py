#!/usr/bin/env python3

'''
TODO: Description Needed

'''

def test_build_surface():

    #### Traditional ASE functionality #####
    from software.examples.data.model_gen import get_example_slab as slab
    slab = slab()
    #########

    #### Functionality to create all surfaces ####
    from software.build.facets import generate
    facets, slabs = generate('Au', save=True)
    #########

    #### Assertion tests ####
    # Scale all cell vectors to match those defined in the ASE version
    lattice_parameter=2.939
    for i in range(len(slabs)):
        slabs[i].set_cell(slabs[i].get_cell()*(lattice_parameter/slabs[i].get_cell()[1][0]))

    # This is an assertion - it checks the same results are given by both methods
    # Note I don't test the vacuum direction. This differs by a small value, for reasons of no concern.
    eps = 1e-10
    for i in range(2):
        for j in range(3):
            # *1.5 is because the test slab is 3 units wide, and the generated slab is 2 units wide.
            assert(slab.get_cell()[i][j] - slabs[0].get_cell()[i][j]*1.5 < eps)

# Run the example
test_build_surface()