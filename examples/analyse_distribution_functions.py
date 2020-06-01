'''
This is an example of usage of the distribution function scripts"

This is useful for looking at the distributions of atoms from a given atom in a system
TODO: Add assertion tests
'''

def test_analyse_distribution_function():
    from carmm.analyse.distribution_functions import distance_distributon_function

    #Build a model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    distance_distributon_function(slab)

test_analyse_distribution_function()
