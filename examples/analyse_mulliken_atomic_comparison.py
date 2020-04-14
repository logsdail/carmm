#!/usr/bin/env python3

import matplotlib.pyplot as plt
from software.analyse.mulliken import parse_mulliken_file
from software.analyse.graphs import get_graph_linetype, get_graph_colour, set_graph_axes_mulliken
from software.analyse.bonds import get_indices_of_elements

output_files = ['data/CO/co_light.log', 'data/Fe/fe_light.log', 'data/Fe-CO/fe-co_light.log']
mulliken_files = ['data/CO/Mulliken.out', 'data/Fe/Mulliken.out', 'data/Fe-CO/Mulliken.out']

for i in range(len(output_files)):
    # Subplots
    ax1 = plt.subplot(len(output_files),1,i+1)

    # Read in atoms information
    output_file = output_files[i]
    from ase.io import read
    atoms = read(output_file)

    # Read in Mulliken data from file
    mulliken_file = mulliken_files[i]
    mulliken_data = parse_mulliken_file(mulliken_file)

    # Collect all the density of states data to plot
    x, all_data = mulliken_data.get_all_plot_data()

    # Collect the indices for each element we are interested in
    fe_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'fe')
    c_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'c')
    o_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'o')

    # Collect the density of states data to plot as a function of atomic label
    x, fe = mulliken_data.get_plot_data(fe_indices, range(mulliken_data.get_nspin()),
                                        range(mulliken_data.get_nkpts()), 'd')
    x, c = mulliken_data.get_plot_data(c_indices, range(mulliken_data.get_nspin()),
                                        range(mulliken_data.get_nkpts()), 'sp')
    x, o = mulliken_data.get_plot_data(o_indices, range(mulliken_data.get_nspin()),
                                        range(mulliken_data.get_nkpts()), 'sp')

    # Put this at the end so it covers everything else and shows the outline of the DOS correctly
    for sp in range(len(all_data)):
        if sp == 0:
            # Here we are filling between lines, rather than plotting a single line.
            # Hence we provide two y-values to plot between.
            plt.fill_between(x, [ 0 * len(x) ], fe[sp], lw=0, facecolor=get_graph_colour(0), label='Fe', interpolate=True)
            plt.fill_between(x, fe[sp], fe[sp]+c[sp], lw=0, facecolor=get_graph_colour(1), label='C', interpolate=True)
            plt.fill_between(x, fe[sp]+c[sp], fe[sp]+c[sp]+o[sp], lw=0, facecolor=get_graph_colour(2), label='O', interpolate=True)
            plt.plot(x, all_data[sp], lw=1, color='black', ls=get_graph_linetype())
        else: # (sp == 1)
            # Note the slightly different definition of the y-axis data here
            # There is a chance that the element doesn't exist in the system, which returns a DOS of zero,
            # and you cannot take the negative of a list of zeros (it transpires)
            # Instead, the approach given just makes each value in the list its negative value separately
            plt.fill_between(x, [0 * len(x)], [-i for i in fe[sp]], lw=0, facecolor=get_graph_colour(0), interpolate=True)
            plt.fill_between(x, [-i for i in fe[sp]], [-i for i in fe[sp]+c[sp]], lw=0, facecolor=get_graph_colour(1), interpolate=True)
            plt.fill_between(x, [-i for i in fe[sp]+c[sp]], [-i for i in fe[sp]+c[sp]+o[sp]], lw=0, facecolor=get_graph_colour(2), interpolate=True)
            plt.plot(x, -(all_data[sp]), lw=1, color='black', ls=get_graph_linetype())

    # Work to rescale axes. Extracts the maximum y-value
    set_graph_axes_mulliken(plt, x, all_data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

    # Add a legend
    plt.legend()

#### Shoe-horn an assertion test in on final HOMO ####
assert(mulliken_data.get_homo() == -4.22285)
########

# Display the graphs
#plt.show()