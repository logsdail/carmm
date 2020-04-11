#!/usr/bin/env python3

import matplotlib.pyplot as plt
from software.analyse.mulliken import parse_mulliken_file, get_graph_linetype, set_graph_axes

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
        plt.plot(x, s[sp], lw=2, color='red', label='s', ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp], lw=2, color='green', label='p', ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp]+d[sp], lw=2, color='blue', label='d', ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp]+d[sp]+f[sp], lw=2, color='black', label='f', ls=get_graph_linetype())
    else: # (sp == 1)
        plt.plot(x, -(s[sp]), lw=2, color='red', ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]), lw=2, color='green', ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]+d[sp]), lw=2, color='blue', ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]+d[sp]+f[sp]), lw=2, color='black', ls=get_graph_linetype())

# Collect all data, which is necessary for managing the axes
all_data = [ s[0] + p[0] +d[0] +f[0] ]
if len(s) > 1:
    all_data.append(s[1] + p[1] +d[1] +f[1])

# Work to rescale axes. Extracts the maximum y-values
set_graph_axes(plt, x, all_data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

# Add a legend
plt.legend()

# Display the graphs
#plt.show()