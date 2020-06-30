'''
This is an example of usage of the distribution function scripts"

This is useful for looking at the distributions of atoms from a given atom in a system
TODO: Add assertion tests
'''

def test_analyse_distribution_function():
    from carmm.analyse.distribution_functions import distance_distribution_function

    #Build a model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    distance_distribution_function(slab, 0.1)

test_analyse_distribution_function()

def test_analyse_radial_distribution_function():
    from carmm.analyse.distribution_functions import radial_distribution_function

    #Build model
    from data.model_gen import get_example_slab as slab_2
    slab_2 = slab_2(adsorbate=True)

    radial_distribution_function(slab_2, 10, 0)


test_analyse_distribution_function()