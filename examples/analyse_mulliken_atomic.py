#!/usr/bin/env python3

import matplotlib.pyplot as plt
from software.analyse.mulliken import parse_mulliken_file, get_graph_linetype, get_indices_of_elements, set_graph_axes

# Read in atoms information
output_file = "data/Fe-CO/fe-co_light.log"
from ase.io import read
atoms = read(output_file)

# Read in Mulliken data from file
mulliken_file = "data/Fe-CO/Mulliken.out"
with open(mulliken_file, 'r') as read_stream:
    lines = read_stream.readlines()

# Parse data from Mulliken file
mulliken_data = parse_mulliken_file(lines)

#### Assertion statements ####
assert(mulliken_data.get_natoms() == 3)
assert(mulliken_data.get_nspin() == 2)
assert(mulliken_data.get_nkpts() == 1)
assert(mulliken_data.get_nstates() == 32)
assert(mulliken_data.get_data_integrity())
#####

# Collect all the density of states data to plot
x, all_data = mulliken_data.get_all_plot_data()

# Collect the indices for each element we are interested in
fe_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'fe')
c_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'c')
o_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'o')

# Collect the density of states data to plot as a function of atomic label
x, fe = mulliken_data.get_plot_data(fe_indices, range(mulliken_data.get_nspin()),
                                    range(mulliken_data.get_nkpts()), 'all')
x, c = mulliken_data.get_plot_data(c_indices, range(mulliken_data.get_nspin()),
                                    range(mulliken_data.get_nkpts()), 'all')
x, o = mulliken_data.get_plot_data(o_indices, range(mulliken_data.get_nspin()),
                                    range(mulliken_data.get_nkpts()), 'all')

# Put this at the end so it covers everything else and shows the outline of the DOS correctly
for sp in range(len(all_data)):
    if sp == 0:
        plt.plot(x, fe[sp], lw=2, color='red', label='Fe', ls=get_graph_linetype())
        plt.plot(x, fe[sp]+c[sp], lw=2, color='green', label='C', ls=get_graph_linetype())
        plt.plot(x, fe[sp]+c[sp]+o[sp], lw=2, color='blue', label='O', ls=get_graph_linetype())
        plt.plot(x, all_data[sp], lw=2, color='black', ls=get_graph_linetype())
    else: # (sp == 1)
        plt.plot(x, -(fe[sp]), lw=2, color='red', ls=get_graph_linetype())
        plt.plot(x, -(fe[sp]+c[sp]), lw=2, color='green', ls=get_graph_linetype())
        plt.plot(x, -(fe[sp]+c[sp]+o[sp]), lw=2, color='blue', ls=get_graph_linetype())
        plt.plot(x, -(all_data[sp]), lw=2, color='black', ls=get_graph_linetype())

# Work to rescale axes. Extracts the maximum y-value
set_graph_axes(plt, x, all_data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

# Add a legend
plt.legend()

# Display the graphs
plt.show()