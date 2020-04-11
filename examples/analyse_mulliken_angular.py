#!/usr/bin/env python3

import matplotlib.pyplot as plt
from software.analyse.mulliken import parse_mulliken_file
from software.analyse.graphs import get_graph_linetype, get_graph_colour, set_graph_axes_mulliken

# Read in data from file
file = "data/CO/Mulliken.out"
with open(file, 'r') as read_stream:
    lines = read_stream.readlines()

# Parse data from Mulliken file
mulliken_data = parse_mulliken_file(lines)

#### Assertion statements ####
assert(mulliken_data.get_natoms() == 2)
assert(mulliken_data.get_nspin() == 2)
assert(mulliken_data.get_nkpts() == 1)
assert(mulliken_data.get_nstates() == 13)
assert(mulliken_data.get_data_integrity())
#####

# Collect the density of states data to plot as a function of angular momenta
x, s = mulliken_data.get_s_plot_data()
x, p = mulliken_data.get_p_plot_data()
x, d = mulliken_data.get_d_plot_data()
x, f = mulliken_data.get_f_plot_data()

# Put this at the end so it covers everything else and shows the outline of the DOS correctly
for sp in range(len(s)):
    if sp == 0:
        plt.plot(x, s[sp], lw=2, label='s', color=get_graph_colour(0), ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp], lw=2, label='p', color=get_graph_colour(1), ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp]+d[sp], lw=2, label='d', color=get_graph_colour(2), ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp]+d[sp]+f[sp], lw=2, label='f', color='black', ls=get_graph_linetype())
    else: # (sp == 1)
        plt.plot(x, -(s[sp]), lw=2, color=get_graph_colour(0), ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]), lw=2, color=get_graph_colour(1), ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]+d[sp]), lw=2, color=get_graph_colour(2), ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]+d[sp]+f[sp]), lw=2, color='black', ls=get_graph_linetype())

# Collect all data, which is necessary for managing the axes
all_data = [ s[0] + p[0] +d[0] +f[0] ]
if len(s) > 1:
    all_data.append(s[1] + p[1] +d[1] +f[1])

# Work to rescale axes. Extracts the maximum y-values
set_graph_axes_mulliken(plt, x, all_data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

# Add a legend
plt.legend()

# Display the graphs
# plt.show()