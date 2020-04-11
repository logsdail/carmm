#!/usr/bin/env python3

import matplotlib.pyplot as plt
from software.analyse.mulliken import parse_mulliken_file, get_graph_linetype

# Read in data from file
file = "data/CO/Mulliken.out"
with open(file, 'r') as read_stream:
    lines = read_stream.readlines()

# Parse data from Mulliken file
mulliken_data = parse_mulliken_file(lines)

# Debug statements to check we have read everything correctly
print("Atoms    : ", mulliken_data.get_natoms())
print("Spins    : ", mulliken_data.get_nspin())
print("K-points : ", mulliken_data.get_nkpts())
print("States   : ", mulliken_data.get_nstates())
print("Data Integrity : ", mulliken_data.get_data_integrity())

#TODO: homo, lumo = atoms.get_homo_lumo_data()
# Collect the data to plot as a function of angular momenta
x, s = mulliken_data.get_s_plot_data()
x, p = mulliken_data.get_p_plot_data()
x, d = mulliken_data.get_d_plot_data()
x, f = mulliken_data.get_f_plot_data()

# Put this at the end so it covers everything else and shows the outline of the DOS correctly
for sp in range(len(s)):
    if sp == 0:
        plt.plot(x, s[sp], lw=2, color='red', label=file, ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp], lw=2, color='green', label=file, ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp]+d[sp], lw=2, color='blue', label=file, ls=get_graph_linetype())
        plt.plot(x, s[sp]+p[sp]+d[sp]+f[sp], lw=2, color='black', label=file, ls=get_graph_linetype())
    else: # (sp == 1)
        plt.plot(x, -(s[sp]), lw=2, color='red', label=file, ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]), lw=2, color='green', label=file, ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]+d[sp]), lw=2, color='blue', label=file, ls=get_graph_linetype())
        plt.plot(x, -(s[sp]+p[sp]+d[sp]+f[sp]), lw=2, color='black', label=file, ls=get_graph_linetype())

# Work to rescale axes. Extracts the maximum y-value
ymax = max(map(max, s+p+d+f))*1.1
ymin = 0
if len(s) > 1:
    ymin = -ymax
    plt.axhline(y=0, xmin=min(x), xmax=max(x), color='black', lw=2)
plt.ylim(ymin, ymax)
plt.yticks([])
plt.ylabel('Density of States (1/eV)')

# Organise x-axis
#plt.xlim(min_point+10, max_point-10)
plt.xlabel(mulliken_data.get_graph_xlabel())

# HOMO
#plt.axvline(x=homo, ymin=-100, ymax=100, color='black', lw=2, ls=line_types[k]) # MFI

# Display the graphs
#plt.show()