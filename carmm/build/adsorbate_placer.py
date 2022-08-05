def rotate_and_place_adsorbate(atoms_ads, atoms_site, bond_length,
                               ads_idx, site_idx, neighb_idx=0, rotation=[0, 0, 0],
                               lps=1, lp_idx=0):
    """

    Place and rotates the adsorbate along a bond length defined by vector -
    vector calculated assuming VSEPR model with one remaining lone pair.
    Rotation takes places along the bonding direction. See documentation for
    rotation variable for directions of rotation.

    Args:
        atoms_ads: Atoms object
            Contains positions of adsorbate.
        atoms_site: Atoms object
            Contains positions of the adsorption site.
        bond_length: float
            Length of the bond for adsorption.
        ads_idx: integer
            Index of the adsorbate atom bonding adsorbate and adsorption site.
        site_idx: integer
            Index of the site atom bonding adsorbate and adsorption site.
        neighb_idx: integer
            Index of neighbour used to form principle axes of rotation.
        rotation: list
            Angles in degrees of rotation for absorbate.
            z - along ads and site bond, x - perpendicular to neighbour and bond, y - perpendicular to x and z.

    Returns:
        ads_and_site: numpy array, (:,3)
            Array of positions of the adsorbed species and site.
        atom_ads: numpy array, (:,3)
            Array of positions for the rotated adsorbed species only.

    """

    # Zeroes and rotates the molecule along the site normal
    junk, zeroed_adsorbate = place_adsorbate(atoms_ads, atoms_site, ads_idx,
                                             site_idx, bond_length, lps=lps, lp_idx=lp_idx)

    # Find appropriate orthogonal rotation axes.
    x_axis, y_axis, z_axis = find_adsorbate_rotation_axes(atoms_site, site_idx, neighb_idx)

    zeroed_adsorbate.rotate(rotation[0], x_axis, center=zeroed_adsorbate.positions[ads_idx])
    zeroed_adsorbate.rotate(rotation[1], y_axis, center=zeroed_adsorbate.positions[ads_idx])
    zeroed_adsorbate.rotate(rotation[2], z_axis, center=zeroed_adsorbate.positions[ads_idx])

    ads_and_site = atoms_site + zeroed_adsorbate

    return ads_and_site, zeroed_adsorbate


def place_adsorbate(atoms_ads, atoms_site, ads_idx, site_idx, bond_length, lps=1, lp_idx=0):
    """
    Place and rotates the adsorbate along a bond length defined by vector.
    Molecule rotated such that the center of positions for the molecule points along
    the bonding direction. Equivalent to rotate_and_place_adsorbate(...,rotation=[0,0,0]).

    Args:
        atoms_ads: Atoms object
            Contains positions of adsorbate.
        atoms_site: Atoms object
            Contains positions of the adsorption site.
        ads_idx: integer
            Index of the adsorbate atom bonding adsorbate and adsorption site.
        site_idx: integer
            Index of the site atom bonding adsorbate and adsorption site.
        bond_length: float
            Length of the bond for adsorption.

    Returns:
        ads_and_site: numpy array, (:,3)
            Array of positions of the adsorbed species and site.
        atom_ads: numpy array, (:,3)
            Array of positions for the rotated adsorbed species only.

    """

    from carmm.analyse.neighbours import neighbours
    import numpy as np

    assert lps < 3,  "Lone pairs greater than 2 not yet implemented."
    assert lps != 0, "No valid adsorption site available."
    assert len(atoms_site) == 0, "Adsorbate site should have at least one other atom attached."

    # Otherwise, just find normal based on average of neighbour vectors.
    # If multiple lps to adsorb to - generate list of potential sites.
    site_normal = find_site_normal(atoms_site, site_idx)

    # Code and logic is sloppy and will be improved.
    if lps>1:

        neighb_list, shell_list = neighbours(atoms_site, [site_idx], 1)

        assert 1 < len(shell_list[1]) < 5, "Site either has too few or too many neighbours VSEPR."
#        assert len(shell_list[1]) > 2, "Not implemented LPs above 2."

        if len(shell_list[1])==2:
            neighb1 = find_generic_normal(atoms_site, site_idx, shell_list[1][0])
            neighb2 = find_generic_normal(atoms_site, site_idx, shell_list[1][1])

            y_rot = np.cross(neighb1, neighb2)
            z_rot = np.cross(y_rot, site_normal)
            print(z_rot)

            if lp_idx==0:
               theta = np.pi/180 * -52
            elif lp_idx==1:
                theta = np.pi/180 * 52

            rot_matrix = normal_rotation_matrix(theta, z_rot)

            site_normal = np.dot(rot_matrix, site_normal.T).T

    if len(atoms_ads.numbers) > 1:
        adsorbate_normal = find_adsorbate_normal(atoms_ads, ads_idx)

        # Returns the zeroed coordinates about atom[ads_idx].
        atoms_ads.positions = rotate_ads2site_vec(atoms_ads, ads_idx, site_normal, adsorbate_normal)
    else:
        atoms_ads.positions = atoms_ads.positions - atoms_ads.positions

    # Places zeroed coordinates at the bonding position
    atoms_ads.positions = atoms_ads.positions + (bond_length * site_normal) + atoms_site.positions[site_idx]

    ads_and_site = atoms_ads + atoms_site

    return ads_and_site, atoms_ads


def find_site_normal(atoms, index):
    """

    Returns a normalised vector specifying the direction of the adsorbate-adsorption
    site bond.

    Args:
        atoms: Atoms object
            Contains the atomic position of the site.
        index: integer
            The atomic index of the desired site atom.

    Returns:
        site_normal: numpy array, (3).
            Contains the normalised direction of the ads+site bond about (0,0,0).

    """
    from carmm.analyse.neighbours import neighbours
    import numpy as np

    neighbour_atoms, shell_list = neighbours(atoms, [index], 1)

    vectors = atoms.positions[neighbour_atoms] - atoms.positions[index]

    site_normal = np.sum(vectors, axis=0)
    site_normal = -site_normal / np.linalg.norm(site_normal)

    return site_normal


def find_adsorbate_normal(atoms, index):
    """

    Returns a normalised vector specifying the direction of the desired bonding
    atom of the adsorbate, and the mean of the cartesian positions of the adsorbate
    molecule.

    Args:
        atoms: Atoms object
            Contains the atomic position of the adsorbate.
        index: integer
            The atomic index of the bonding atom.

    Returns:
        site_normal: numpy array, (3).
            Contains the normalised direction of bonding atom->center of positions
            about (0,0,0).

    """

    import numpy as np

    ads_com = atoms.get_positions().mean(axis=0)

    molecule_vector = (atoms.positions[index] - ads_com) / np.linalg.norm(atoms.positions[index] - ads_com)

    return molecule_vector

def find_generic_normal(atoms, index1, index2):
    """
    Returns a normalised vector specifying the direction of a generic
    bond.

    Args:
        atoms: Atoms object
            Contains the atomic position of the adsorbate.
        index1: integer
            The atomic index of atom 1.
        index2: integer
            The atomic index of atom 2.

    Returns:
        site_normal: numpy array, (3).
            Contains the normalised direction of atom1 to atom2.
            about (0,0,0).

    """

    import numpy as np

    vector = atoms.positions[index2] - atoms.positions[index1]

    molecule_vector = (vector) / np.linalg.norm(vector)

    return molecule_vector

def find_adsorbate_rotation_axes(atoms_site, site_idx, neighb_idx):
    """

    Find the x, y and z rotation axes, defined by the direction of the bond
    to the adsorbate (site_normal) and the direction of the bond to a given
    pre-existing neighbour atom. If only one around the adsorbate site is found,
    the next atom its neighbours are used.

    Args:
        atoms_site: Atoms object
            Contains positions of the adsorption site.
        site_idx: integer
            Index of the site atom bonding adsorbate and adsorption site.
        neighb_idx: integer
            Index of neighbour used to form principle axes of rotation.

    Returns:
        x_axis: numpy array, (3)
            x rotation axis (perpendicular to site bonds to the adsorbate and
            another neighbour atom).
        y_axis, numpy array, (3)
            y rotation axis (perpendicular given x and y rotation axes).
        z_axis, numpy array, (3)
            z rotation axis (defined by the site normal).

    """
    import numpy as np
    from carmm.analyse.neighbours import neighbours

    # Rotations performed on a RHS axis, with rotations about site normal z, x out plane wrt.
    # bond vector of another neighbour to the shell site and z, and y perpendicular to x and z.
    site_norm = find_site_normal(atoms_site, site_idx)
    neighbour_atoms, shell_list = neighbours(atoms_site, [site_idx], 1)

    # Clean up list for clarity.
    neighbour_atoms.remove(site_idx)
    assert len(neighbour_atoms) > 0, "Adsorption site does not have any neighbours. Cannot define rotation axes."

    # If the site only has one neighbour, the neighbour normal is unobtainable.
    # Next site is used to find an appropriate x and y rotation axis.
    if len(neighbour_atoms) == 1:
        new_site_idx = neighbour_atoms[0]
        neighbour_atoms, shell_list = neighbours(atoms_site, [new_site_idx], 1)

        # Clean up list - removes error producing options.
        neighbour_atoms.remove(new_site_idx)
        neighbour_atoms.remove(site_idx)

        # Order list for compatability reasons.
        neighbour_atoms.sort()

        vectors = atoms_site.positions[neighbour_atoms] - atoms_site.positions[new_site_idx]
        neighb_vector = vectors[neighb_idx] / np.linalg.norm(vectors[neighb_idx])

    else:
        # Order list for compatability reasons.
        neighbour_atoms.sort()

        vectors = atoms_site.positions[neighbour_atoms] - atoms_site.positions[site_idx]
        neighb_vector = vectors[neighb_idx] / np.linalg.norm(vectors[neighb_idx])

    # Given the choice of two orthogonal vectors wrt the site normal is arbitrary, we choose the cross product
    # with respect to the vector from site to its neighbours, then the cross product with this newly generated
    # axis with the site normal.
    z_axis = site_norm
    x_axis = np.cross(neighb_vector, site_norm)
    y_axis = np.cross(x_axis, z_axis)

    return x_axis, y_axis, z_axis


def rotate_ads2site_vec(atoms, atoms_idx, site_vec, mol_vec):
    """

    Performs a rotation of the molecule onto the site normal.
    ROTATION MATRIX TAKEN FROM https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724

    Args:
        atoms: Atoms object
            Contains positions of adsorbate.
        atoms_idx: integer
            Index of the adsorbate atom bonding adsorbate and adsorption site.
        site_vec: numpy array, size=(3)
            vector defining the direction
        mol_vec: numpy array, size=(3)
            vector defining direction from the bonding atom to the mean of
            the atomic positions.

    Returns:
        rotated_positions: numpy array, size=(Nx3)
            Position of the atoms rotated so site_vec and mol_vec align.

    """

    import numpy as np

    theta = np.arccos(np.dot(mol_vec, site_vec))

    axis = np.cross(mol_vec, site_vec)

    axis = axis / np.linalg.norm(axis)

    rotation_matrix = -normal_rotation_matrix(theta, axis)

    zeroed_atomco = atoms.positions - atoms.positions[atoms_idx]
    rotated_positions = np.dot(rotation_matrix, zeroed_atomco.T).T

    return rotated_positions

def normal_rotation_matrix(theta, ax):
    """

    Returns a (3x3) rotation matrix for a rotation about a given
    axis by a theta value.

    Args:
        theta: float
            Angle of rotation in radians.
        ax: numpy array, size = (3)
            Direction of the rotation axis, centered on [0, 0, 0]

    Returns:
        rotation_matrix: numpy array, size = (3x3)
            Rotation matrix to be applied to atomic positions.

    """

    import numpy as np

    sina = np.sin(theta)
    cosa = np.cos(theta)

    onemincosa = 1.0 - cosa

    rotation_matrix = np.array([[(ax[0] * ax[0] * onemincosa) + cosa,
                                 (ax[1] * ax[0] * onemincosa) - (sina * ax[2]),
                                 (ax[2] * ax[0] * onemincosa) + (sina * ax[1])],
                                [(ax[0] * ax[1] * onemincosa) + (sina * ax[2]),
                                 (ax[1] * ax[1] * onemincosa) + cosa,
                                 (ax[2] * ax[1] * onemincosa) - (sina * ax[0])],
                                [(ax[0] * ax[2] * onemincosa) - (sina * ax[1]),
                                 (ax[1] * ax[2] * onemincosa) + (sina * ax[0]),
                                 (ax[2] * ax[2] * onemincosa) + cosa]])

    return rotation_matrix
