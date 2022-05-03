#!/usr/bin/env python3

'''

An example showing how to plot the total volume unoccupied by molecular vdW spheres.

A WORD OF WARNING - this script shows the total volume that can be occupied by a probe sphere
of a set radius, but does not differentiate between volumes that are accessible/inaccessible
to the probe. More sophisticated techniques such as Delauney triangulation required.

'''

def test_void_find():

    from carmm.analyse.meshgrid.meshgrid_meshobject import Mesh
    from carmm.analyse.meshgrid.meshgrid_functions import atom_mesh_build_mask, distance_meshgrid2point
    from carmm.analyse.meshgrid.meshgrid_void import void_find, void_build_mask, void_analysis
    import matplotlib.pyplot as plt
    import numpy as np
    from ase.io import read

    # Initialise the Mesh object, with a default underlying meshgrid of 50x50x50 grid points.
    mfi_mesh = Mesh(np.array([20.22614449, 19.82125040, 13.36948553, 90, 90, 90]),
                    nx = 20, ny = 20, nz = 20, pbc=[1, 1, 1])

    # Read atomic coordinates and create a meshgrid using periodic boundary conditions
    # Meshgrid arrays print the x, y, z coordinates for all given points into three (nx,ny,nz)
    # arrays (xx, yy and zz respectively).
    atom = read('data/Zeolite/MFI_framework_geom.xyz')
    mol_xx, mol_yy, mol_zz = atom_mesh_build_mask(mfi_mesh, atom)

    # First find the maximum volume sphere not containing an atom around each grid (probe) point.
    # Coarseness controls the number of probe points by increasing the interval between each probe point.
    # Outputs a list of the void center coordinates and their maximum radii.
    void_centers, void_radii = void_find(mfi_mesh, atom, coarseness=2)

    # Plot the void centers onto the grid. Min_void prevents spheres of radii smaller than a given
    # value being plotted.
    void_xx,void_yy,void_zz=void_build_mask(mfi_mesh, void_centers, void_radii, min_void=1.4)

    # Print basic output details such as the total void volume and the maximum void volume,
    # radius and position.
    # Only requires one void meshgrid to determine whether a site is unoccupied/occupied.
    periodic_volume = void_analysis(mfi_mesh, void_centers, void_radii, void_xx)

    mfi_mesh_obc = Mesh(np.array([20.22614449, 19.82125040, 13.36948553, 90, 90, 90]),
                        nx = 20, ny = 20, nz = 20, pbc=[0, 0, 0])

    test_x, test_y, test_z = atom.positions[0][0], atom.positions[0][1], atom.positions[0][2] 
    periodic_dist = distance_meshgrid2point(test_x, test_y, test_z, mfi_mesh)
    obc_dist      = distance_meshgrid2point(test_x, test_y, test_z, mfi_mesh_obc)

    print(np.mean(periodic_dist))
    print(np.mean(obc_dist))
    
    assert(np.mean(periodic_dist) < np.mean(obc_dist))
    assert(np.isclose(periodic_volume, 890.4178197211703, atol=1e-4, rtol=0.0))

test_void_find()
