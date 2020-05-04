#!/usr/bin/env python3

'''
This is a simple example and QA test for extracting and plotting Mulliken data from FHI-aims,
specifically for a bulk material with k-points in the Mulliken data
Note that this routine plots a graph with only one spin channel, as the calculation is spin-paired

This is useful when trying to understand the electronic structure of your model
'''

def test_analyse_mulliken_atomic_angular():

    import matplotlib.pyplot as plt
    from carmm.analyse.mulliken import parse_mulliken_file
    from carmm.analyse.graphs import get_graph_linetype, get_graph_colour, set_graph_axes_mulliken
    from carmm.analyse.bonds import get_indices_of_elements

    # Read in atoms information
    output_file = "data/TiO2/tio2_rutile_light.log"
    from ase.io import read
    atoms = read(output_file)

    # Read in Mulliken data from file
    mulliken_file = "data/TiO2/Mulliken.out"
    mulliken_data = parse_mulliken_file(mulliken_file)

    #### Assertion statements ####
    assert(mulliken_data.get_natoms() == 6)
    assert(mulliken_data.get_nspin() == 1)
    assert(mulliken_data.get_nkpts() == 184)
    assert(mulliken_data.get_nstates() == 61)
    assert(mulliken_data.get_data_integrity())
    #####

    # Collect all the density of states data to plot
    x, all_data = mulliken_data.get_all_plot_data()

    # Collect the indices for each element we are interested in
    ti_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'ti')
    o_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'o')

    # Collect the density of states data to plot as a function of atomic label
    x, ti_d = mulliken_data.get_plot_data(ti_indices, range(mulliken_data.get_nspin()),
                                        range(mulliken_data.get_nkpts()), 'd')
    x, o_p = mulliken_data.get_plot_data(o_indices, range(mulliken_data.get_nspin()),
                                        range(mulliken_data.get_nkpts()), 'p')

    # Plot the specific cases we are interested in; notice that this is spin-paired
    plt.plot(x, ti_d[0], lw=2, label='Ti 3d', color=get_graph_colour(0), ls=get_graph_linetype())
    plt.plot(x, ti_d[0]+o_p[0], lw=2, label='O 2p', color=get_graph_colour(1), ls=get_graph_linetype())
    # Put this at the end so it covers everything else and shows the outline of the DOS correctly
    plt.plot(x, all_data[0], lw=2, label='All', color='black', ls=get_graph_linetype())

    # Work to rescale axes. Extracts the maximum y-value
    set_graph_axes_mulliken(plt, x, all_data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

    # Add a legend
    plt.legend()

    # Display the graphs
    # plt.show()

# Run the example
test_analyse_mulliken_atomic_angular()