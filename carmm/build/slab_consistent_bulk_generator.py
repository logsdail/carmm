
# This functionality will help to generate a bulk model from a slab.

from ase import Atoms


def is_close_to_integer(arr, tolerance=1e-5):
    """
    Check if elements in the array are close to integers within a given tolerance. This will be used to check if the x-
    and y-coordinate of an atom is a linear combination the 'a' and 'b' cell vector of the slab model

    Args:
    - arr: Numpy array
    - tolerance: Tolerance for closeness to integers (default is 1e-5)

    Returns:
    - Boolean array indicating if elements are close to integers
    """
    import numpy as np

    diff = np.abs(arr - np.round(arr))
    close_to_integer = diff < tolerance
    return close_to_integer


def bulk_identifier(slab: Atoms, cutoff_distance=10.0):
    """
    This function helps to identify the repeating unit in the z-direction which ultimately represents the
    bulk atoms in a slab structure. This new bulk structure will be consistent with the slab model and the
    computed bulk energy for this structure will help to obtain surface energies that are convergent with
    increasing slab thickness. Further detailed information is available in the following literature article.
    Schultz, Peter A.
    "First-principles calculations of metal surfaces. I. Slab-consistent bulk reference for convergent surface
    properties." Physical Review B 103.19 (2021): 195426.

    How this functionality works:
    1. the difference between the x,y and z coordinates of an atom and its corresponding periodic image will be a linear
    combination of the cell vectors a,b and c.
    2. This idea is extended here in case of a slab model where the linear combination check
    is considered only w.r.t to the a and b vectors

    Args:
    - slab: ASE Atoms object representing the slab structure
    - cutoff_distance: Cutoff distance for identifying neighboring atoms (default is 10.0). This will be used to
     perform further checks on atoms which satisfy the condition of linear combination.

    Returns:
    - ASE Atoms object representing the new bulk structure
    """
    from ase.build.tools import sort
    import numpy as np

    # Sort the slab based on z-coordinate
    slab = sort(slab, np.array(slab.get_positions()[:, 2]).tolist())
    a, b, c = slab.cell  # Extract cell vectors
    slab_copy = slab.copy() # make a copy of the slab

    # Calculate all pairwise distances between atoms
    distances = slab.get_all_distances()
    z_coord_values = []
    cutoff_distance = cutoff_distance

    for atom_i in slab:
        index_i = atom_i.index

        for atom_j in slab_copy:
            if atom_i.symbol == atom_j.symbol and atom_i.index != atom_j.index:
                diff_xy = np.array(atom_j.position[:2] - atom_i.position[:2])
                cell_vec_xy = np.array(slab.cell[:2, :2])
                int_vec = np.linalg.solve(cell_vec_xy.transpose(), diff_xy)

                if is_close_to_integer(int_vec).all():
                    # further check is done to see if the coordination environment of the atom j is similar ot atom i
                    # based on cutoff distance.
                    index_j = atom_j.index
                    cutoff_obeyed_i = [dist for ind, dist in enumerate(distances[index_i].tolist())
                                       if
                                       dist <= cutoff_distance and slab[ind].position[2] / atom_i.position[2] >= 0.9999]
                    cutoff_obeyed_j = [dist for ind, dist in enumerate(distances[index_j].tolist())
                                       if
                                       dist <= cutoff_distance and slab[ind].position[2] / atom_j.position[2] >= 0.9999]

                    if len(cutoff_obeyed_j) == len(cutoff_obeyed_i):
                        dist_diff = np.linalg.norm(np.array(cutoff_obeyed_j)) / np.linalg.norm(
                            np.array(cutoff_obeyed_i))
                        if 0.99 <= dist_diff <= 1.01:
                            z_dist = atom_j.position[2] - atom_i.position[2]
                            z_coord_values.append(z_dist)

    z_coord_values = np.array(z_coord_values)
    slab.center()
    # set the cell parameters to the original 'a' and 'b' whereas the 'c' vector is changed to a new value which
    # represents the minimum repeating unit in the z-direction
    slab.set_cell(np.array([a, b, [0, 0, min(np.abs(z_coord_values))]]))

    # Modify positions to lie within the cell
    x_coord = slab.get_positions()[:, 0]
    y_coord = slab.get_positions()[:, 1]
    new_z_coord = slab.get_positions()[:, 2] - min(slab.get_positions()[:, 2])
    slab.set_positions(np.array([x_coord, y_coord, new_z_coord]).transpose())

    # Remove atoms beyond cell boundaries
    atoms_to_del = [atom_i.index for atom_i in slab if atom_i.position[2] > slab.cell[2, 2]]
    new_bulk = slab[[i for i in range(len(slab)) if i not in atoms_to_del]]

    return new_bulk # return the new bulk model