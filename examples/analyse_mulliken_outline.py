#!/usr/bin/env python3

'''
This is a simple example and QA test for extracting and plotting Mulliken data from FHI-aims
Note that this routine plots a graph with both spin-up and spin-down data.

This is useful when trying to understand the electronic structure of your model
'''

def test_mulliken_outline():

    import matplotlib.pyplot as plt
    from carmm.analyse.mulliken import parse_mulliken_file, write_dos_to_csv
    from carmm.analyse.graphs import get_graph_linetype, get_graph_colour, set_graph_axes_mulliken

    # Read in data from file
    file = "data/Fe/Mulliken.out"
    mulliken_data = parse_mulliken_file(file)

    #### Assertion statements ####
    assert(mulliken_data.get_natoms() == 1)
    assert(mulliken_data.get_nspin() == 2)
    assert(mulliken_data.get_nkpts() == 1)
    assert(mulliken_data.get_nstates() == 19)
    assert(mulliken_data.get_data_integrity())
    #####

    # Collect all the density of states data to plot
    x, data = mulliken_data.get_all_plot_data()
    # Get the angular contributions so we can see these too
    x, s = mulliken_data.get_orbital_plot_data(orbital='s')
    x, p = mulliken_data.get_orbital_plot_data(orbital='p')
    x, d = mulliken_data.get_orbital_plot_data(orbital='d')
    x, f = mulliken_data.get_orbital_plot_data(orbital='f')

    # The returned array has either one entries if spin-collinear, or two entries for spin up and down.
    for spin in range(mulliken_data.get_nspin()):
        if spin == 0:
            # Plot the total outline last
            plt.plot(x, data[spin], lw=2, color='black', ls=get_graph_linetype())
        else: # (spin == 1)
            # Plot all the angular contributions to show what is possible
            plt.fill_between(x, [0 * len(s[spin])], -s[spin], lw=2, color=get_graph_colour(0), ls=get_graph_linetype(), label='s')
            plt.fill_between(x, -s[spin], -(s[spin]+p[spin]), lw=2, color=get_graph_colour(1), ls=get_graph_linetype(), label='p')
            plt.fill_between(x, -(s[spin]+p[spin]),-(s[spin]+p[spin]+d[spin]), lw=2, color=get_graph_colour(2), ls=get_graph_linetype(), label='d')
            plt.fill_between(x, -(s[spin]+p[spin]+d[spin]), -(s[spin]+p[spin]+d[spin]+f[spin]), lw=2, color=get_graph_colour(3), ls=get_graph_linetype(), label='f')
            # Plot the total outline last
            plt.plot(x, -(data[spin]), lw=2, color='black', ls=get_graph_linetype())

    # Work to rescale axes. Extracts the maximum y-value
    set_graph_axes_mulliken(plt, x, data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

    # Example of how to save data to csv file for exporting
    write_dos_to_csv('all_data.csv', x, data)

    # Display the graphs
    plt.legend()
    #plt.show()

# Run the example/test
test_mulliken_outline()