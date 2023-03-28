'''
This is an example of usage of the distribution function scripts

This is useful for looking at the distributions of atoms from a given atom in a system
'''

def test_analyse_radial_distribution_function():
    from carmm.analyse.distribution_functions import radial_distribution_function

    #Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    distances = radial_distribution_function(model=slab, radius=10, position=0, verbose=True)
    assert(len(distances) == 88)
    assert(1e-5 > abs(distances[0] - 2.938999))
    assert(1e-5 > abs(distances[87] - 9.936709))

    from carmm.analyse.distribution_functions import plot_distribution_function
    plt = plot_distribution_function(distances, title='Radial Distribution Function')
    #plt.show()

def test_analyse_element_radial_distribution_function():
    from carmm.analyse.distribution_functions import element_radial_distribution_function
    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    #Get chemical symbols
    chemical_symbols = slab.get_chemical_symbols()

    distances = element_radial_distribution_function(model=slab, radius=10, element=chemical_symbols[0], verbose=True)
    assert (len(distances) == 78)
    assert (1e-5 > abs(distances[0] - 2.938999))
    assert (1e-5 > abs(distances[77] - 9.74756))

    from carmm.analyse.distribution_functions import plot_distribution_function
    plt = plot_distribution_function(distances, title=chemical_symbols[0] + ' Radial Distribution Function')
    # plt.show()

def test_analyse_average_distribution_function():
    from carmm.analyse.distribution_functions import average_distribution_function

    #Build model
    from data.model_gen import get_example_slab as slab
    slab_1 = slab(adsorbate=True)
    slab_2 = slab(adsorbate=True)
    slab_2.positions *= 1.05
    slab_trajectory = [slab_1, slab_2]

    all_data, snapshots = average_distribution_function(trajectory=slab_trajectory, samples=2)
    assert(len(all_data) == 420)
    assert(1e-5 > abs(sum(all_data)-1754.778323))

    # Example how to plot data
    from carmm.analyse.distribution_functions import plot_distribution_function
    from matplotlib import pyplot as plt
    # Clean axes before new plot
    plt.cla()

    # Plot all snapshots
    for i in range(len(snapshots)):
        plot_distribution_function(snapshots[i], color='#808080', density=True)
    # Plot mean (i.e. average bins over all data, rather than separate snapshots).
    # This probably need double checking with something known e.g. H2O TODO
    plot_distribution_function(all_data, title='Last ' + str(len(snapshots)) + ' Snapshots', density=True,
                               # Additional graph plotting variables
                               label='Mean', color="#000000")
    plt.legend()
    #plt.show()

def test_radius_of_gyration():
    from carmm.analyse.distribution_functions import radius_of_gyration

    # Build model
    from ase.build import molecule
    water = molecule('H2O')

    # Get radius of gyration
    rog = radius_of_gyration(model=water)
    assert(1e-5 > abs(rog-0.317063))




def test_rdf():
    from carmm.analyse.distribution_functions import rdf

    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    element1 = 'Au'
    element2 = 'Au'
    rdf = rdf(slab, element1, element2)

    from carmm.analyse.distribution_functions import plot_rdf

    y_lab = '$g(r)_{Au-Au}$'
    plt = plot_rdf(rdf, y_lab=y_lab)
    plt.show()


# test_analyse_radial_distribution_function()
# test_analyse_element_radial_distribution_function()
# test_analyse_average_distribution_function()
# test_radius_of_gyration()

test_rdf()