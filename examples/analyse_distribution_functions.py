'''
This is an example of usage of the distribution function scripts"

This is useful for looking at the distributions of atoms from a given atom in a system
TODO: Add assertion tests
       - need to add an example for average_distribution_function - will require creating an MD traj
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
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    radial_distribution_function(slab, 10, 0)

def test_analyse_average_distribution_function():
    from carmm.analyse.distribution_functions import average_distribution_function

    #Build model
    from data.model_gen import get_example_slab as slab
    slab_1 = slab(adsorbate=True)
    slab_2 = slab(adsorbate=True)
    slab_2.positions *= 1.05
    slab_trajectory = [slab_1, slab_2]

    average_distribution_function(slab_trajectory, 2)

#test_analyse_distribution_function()
#test_analyse_radial_distribution_function()
test_analyse_average_distribution_function()