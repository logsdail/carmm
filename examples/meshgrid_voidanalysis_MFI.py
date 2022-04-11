#!/usr/bin/env python3

from carmm.meshgrid.meshgrid_unitcell import UnitCell
from carmm.meshgrid.meshgrid_functions import atom_mesh_build_mask
from carmm.meshgrid.meshgrid_void import void_find
from carmm.meshgrid.meshgrid_void import void_build_mask
from carmm.meshgrid.meshgrid_void import void_analysis
import matplotlib.pyplot as plt
from ase import io

# An example showing how to plot the total volume not occupied by molecular vdW volume.

# A WORD OF WARNING - this script shows the total volume that can be occupied by a probe sphere
# of a set radius, but does not differentiate between volumes that are accessible/inaccessible
# to the probe. More sophisticated techniques such as Delauney triangulation required.

# Initialise the unit cell object, with a default underlying meshgrid of 50x50x50 grid points.
ucell=UnitCell()
ucell.define_unit_cell(20.22614449,19.82125040,13.36948553)

# Read atomic coordinates and create a meshgrid using periodic boundary conditions
# (controlled by minimum image convention (MIC) true/false).
# Meshgrid arrays print the x, y, z coordinates for all given points into three (nx,ny,nz)
# arrays (mol_xx, mol_yy and mol_zz respectively).
atom=io.read('data/Zeolite/MFI_framework_geom.xyz')
mol_xx,mol_yy,mol_zz=atom_mesh_build_mask(ucell,atom,mic=True)

# First find the maximum volume sphere not containing an atom around each grid (probe) point.
# Coarseness controls the number of probe points by increasing the interval between each probe point.
# Outputs a list of the void center coordinates and their maximum radii.
void_centers,void_radii=void_find(ucell, atom, mic=True, coarseness=2)

# Plot the void centers onto the grid. Min_void prevents spheres of radii smaller than a given
# value being plotted.
void_xx,void_yy,void_zz=void_build_mask(ucell, void_centers, void_radii, mic=True, min_void=1.4)

# Print basic output details such as the total void volume and the maximum void volume, radius and position.
# Only requires one void meshgrid to determine whether a site is unoccupied/occupied.
void_analysis(ucell, void_centers, void_radii, void_xx)

# Plots the meshgrid for the molecules and void.
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(mol_xx,mol_yy,mol_zz)
ax.scatter(void_xx,void_yy,void_zz)

ax.set_xlim(0,ucell.dim[0])
ax.set_ylim(0,ucell.dim[1])
ax.set_zlim(0,ucell.dim[2])
plt.show()