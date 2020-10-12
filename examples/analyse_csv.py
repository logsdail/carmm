#!/usr/bin/env python3

'''
This example shows how to load datasets from CSV format whereby the a heatmap is then to be plotted.

This comes in useful when working with complex phase diagrams
'''

def test_analyse_csv():
    from carmm.analyse.graphs import load_xyz_data_from_csv, set_graph_axes_heatmap
    import matplotlib.pyplot as plt

    # Load and triangulate the data
    x, y, z = load_xyz_data_from_csv('data/csv/data_scan.dat')

    # Plot the data
    plt.tricontourf(x, y, z, 1000)
    set_graph_axes_heatmap(plt, x, y)

    #### ASSERTION TEST ####
    assert(plt.xlim()[0] == 6.1 and plt.xlim()[1] == 6.26)
    assert(plt.ylim()[0] == 90.0 and plt.ylim()[1] == 95.13)
    #########

    # Visualise plot
    #plt.show()

# Run example
test_analyse_csv()