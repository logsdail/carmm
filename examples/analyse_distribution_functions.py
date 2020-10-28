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

    distances = distance_distribution_function(slab, 0.1, plot=True)
    # Check values are stil lthe same and ordering is correct
    assert(len(distances) == 210)
    assert(1e-5 > abs(distances[0] - 1.178657))
    assert(1e-5 > abs(distances[-1] - 7.132005))

def test_analyse_radial_distribution_function():
    from carmm.analyse.distribution_functions import radial_distribution_function

    #Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    distances = radial_distribution_function(slab, 10, 0, plot=True)
    # There is an error here. The supercell isn't big enough as we should have more than 22 interactions!
    # TODO: Make the sphere cutout procedure more robust i.e. check size of cell against radius
    assert(len(distances) == 22)
    assert(1e-5 > abs(distances[0] - 0.000000))
    assert(1e-5 > abs(distances[-1] - 9.747560))

def test_analyse_average_distribution_function():
    from carmm.analyse.distribution_functions import average_distribution_function

    #Build model
    from data.model_gen import get_example_slab as slab
    slab_1 = slab(adsorbate=True)
    slab_2 = slab(adsorbate=True)
    slab_2.positions *= 1.05
    slab_trajectory = [slab_1, slab_2]

    average_distribution_function(slab_trajectory, 2, plot=True)

test_analyse_distribution_function()
test_analyse_radial_distribution_function()
test_analyse_average_distribution_function()