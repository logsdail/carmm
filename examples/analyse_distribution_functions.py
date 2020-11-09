'''
This is an example of usage of the distribution function scripts"

This is useful for looking at the distributions of atoms from a given atom in a system
TODO: Add assertion tests
       - need to add an example for average_distribution_function - will require creating an MD traj
'''

def test_analyse_distribution_function():
    from carmm.analyse.distribution_functions import radial_distribution_function

    #Build a model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    distances = radial_distribution_function(slab)
    # Check values are stil lthe same and ordering is correct
    assert(len(distances) == 210)
    assert(1e-5 > abs(distances[0] - 1.178657))
    assert(1e-5 > abs(distances[-1] - 7.132005))

    # Test plot
    from carmm.analyse.distribution_functions import plot_distribution_function
    plt = plot_distribution_function(distances, title='Radial Distribution Function')
    #plt.show()

def test_analyse_radial_distribution_function():
    from carmm.analyse.distribution_functions import extended_radial_distribution_function

    #Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    distances = extended_radial_distribution_function(slab, 10, 0)
    assert(len(distances) == 79)
    assert(1e-5 > abs(distances[0] - 2.938999))
    assert(1e-5 > abs(distances[-1] - 9.936709))

    from carmm.analyse.distribution_functions import plot_distribution_function
    plt = plot_distribution_function(distances, title='Radial Distribution Function')
    #plt.show()

def test_analyse_average_distribution_function():
    from carmm.analyse.distribution_functions import average_distribution_function

    #Build model
    from data.model_gen import get_example_slab as slab
    slab_1 = slab(adsorbate=True)
    slab_2 = slab(adsorbate=True)
    slab_2.positions *= 1.05
    slab_trajectory = [slab_1, slab_2]

    mean = average_distribution_function(slab_trajectory, 2, plot=True)
    assert(len(mean) == 420)
    assert(1e-5 > abs(sum(mean)-1754.778323))

test_analyse_distribution_function()
test_analyse_radial_distribution_function()
test_analyse_average_distribution_function()