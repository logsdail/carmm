#!/usr/bin/env python3

import matplotlib.pyplot as plt
from software.analyse.mulliken import parse_mulliken_file, write_to_csv
from software.analyse.graphs import get_graph_linetype, set_graph_axes_mulliken

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

# Put this at the end so it covers everything else and shows the outline of the DOS correctly
for sp in range(len(data)):
    if sp == 0:
        plt.plot(x, data[sp], lw=2, color='black', ls=get_graph_linetype())
    else: # (sp == 1)
        plt.plot(x, -(data[sp]), lw=2, color='black', ls=get_graph_linetype())

# Work to rescale axes. Extracts the maximum y-value
set_graph_axes_mulliken(plt, x, data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

# Example of how to save data to csv file for exporting
# write_to_csv('all_data.csv', x, data)

# Display the graphs
# plt.show()