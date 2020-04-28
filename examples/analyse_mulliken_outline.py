#!/usr/bin/env python3

'''
TODO: Description Needed

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
    x, s = mulliken_data.get_s_plot_data()
    x, p = mulliken_data.get_p_plot_data()
    x, d = mulliken_data.get_d_plot_data()
    x, f = mulliken_data.get_f_plot_data()

    # Put this at the end so it covers everything else and shows the outline of the DOS correctly
    for sp in range(len(data)):
        if sp == 0:
            # Plot the total outline last
            plt.plot(x, data[sp], lw=2, color='black', ls=get_graph_linetype())
        else: # (sp == 1)
            # Plot all the angular contributions to show what is possible
            plt.fill_between(x, [0 * len(s[sp])], -s[sp], lw=2, color=get_graph_colour(0), ls=get_graph_linetype(), label='s')
            plt.fill_between(x, -s[sp], -(s[sp]+p[sp]), lw=2, color=get_graph_colour(1), ls=get_graph_linetype(), label='p')
            plt.fill_between(x, -(s[sp]+p[sp]),-(s[sp]+p[sp]+d[sp]), lw=2, color=get_graph_colour(2), ls=get_graph_linetype(), label='d')
            plt.fill_between(x, -(s[sp]+p[sp]+d[sp]), -(s[sp]+p[sp]+d[sp]+f[sp]), lw=2, color=get_graph_colour(3), ls=get_graph_linetype(), label='f')
            # Plot the total outline last
            plt.plot(x, -(data[sp]), lw=2, color='black', ls=get_graph_linetype())

    # Work to rescale axes. Extracts the maximum y-value
    set_graph_axes_mulliken(plt, x, data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

    # Example of how to save data to csv file for exporting
    write_dos_to_csv('all_data.csv', x, data)

    # Display the graphs
    plt.legend()
    #plt.show()

# Run the example/test
test_mulliken_outline()