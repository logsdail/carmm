def place_adsorbate(atoms_site, atoms_ads, ads_idx, site_idx, bond_length):

    site_normal = find_site_normal(atoms_site, site_idx)

    if len(atoms_ads.numbers)>1:
        adsorbate_normal = find_adsorbate_normal(atoms_ads, ads_idx)

        # Returns the zeroed coordinates about atom[ads_idx].
        atoms_ads.positions = rotate_mol2site_vec(atoms_ads, ads_idx, site_normal, adsorbate_normal)
    else:
        atoms_ads.positions = atoms_ads.positions - atoms_ads.positions

    # Places zeroed coordinates at the bonding position
    atoms_ads.positions = atoms_ads.positions + (bond_length * site_normal) + atoms_site.positions[site_idx]

    ads_and_site = atoms_ads + atoms_site

    return ads_and_site, atoms_ads


def rotate_adsorbate_site_normal(atoms_mol, atoms_site, site_idx, ads_idx, neighb_idx, rotation=[0,0,0]):

    import numpy as np
    from carmm.analyse.neighbours import neighbours

    # Basic wrapper for rotation about site normal.
    # Rotations performed on a RHS axis, with rotations about site normal x, out of plane z, and orthogonal y.

    site_norm = find_site_normal(atoms_site, site_idx)

    neighbour_atoms, shell_list = neighbours(atoms_site, [site_idx], 1)

    # If the site only has one neighbour, the normal is unobtainable.
    # Next site is used to fine an appropriate site.

    if len(neighbour_atoms)<3:
        new_site_idx=neighbour_atoms[1]
        neighbour_atoms, shell_list = neighbours(atoms_site, [new_site_idx], 1)

        vectors=atoms_site.positions[neighbour_atoms]-atoms_site.positions[new_site_idx]
        neighb_vector=vectors[neighb_idx+1]/np.linalg.norm(vectors[neighb_idx+1])

    else:
        vectors=atoms_site.positions[neighbour_atoms]-atoms_site.positions[site_idx]
        neighb_vector=vectors[neighb_idx+1]/np.linalg.norm(vectors[neighb_idx+1])

    # Given the choice of two orthogonal vectors wrt the site normal is arbitrary, we choose the cross product
    # with respect to the vector from site to its neighbours, then the cross product with this newly generated
    # axis with the site normal.
    z_axis=site_norm
    x_axis=np.cross(neighb_vector,site_norm)
    y_axis=np.cross(x_axis,z_axis)

    print(f"Axes: {x_axis,y_axis,z_axis}")
    atoms_mol.rotate(rotation[0],x_axis,center=atoms_mol.positions[ads_idx])
    atoms_mol.rotate(rotation[1],y_axis,center=atoms_mol.positions[ads_idx])
    atoms_mol.rotate(rotation[2],z_axis,center=atoms_mol.positions[ads_idx])

    ads_and_site = atoms_mol + atoms_site

    return ads_and_site, atoms_mol

def find_site_normal(atoms, index):

    from carmm.analyse.neighbours import neighbours
    import numpy as np

    neighbour_atoms, shell_list = neighbours(atoms, [index], 1)

    vectors=atoms.positions[neighbour_atoms]-atoms.positions[index]

    site_normal=np.sum(vectors,axis=0)
    site_normal=-site_normal/np.linalg.norm(site_normal)

    print(site_normal+atoms.positions[index])

    return site_normal

def find_adsorbate_normal(atoms, index):

    import numpy as np

    ads_com = atoms.get_center_of_mass()
    ads_com = atoms.get_positions().mean(axis=0)

    molecule_vector = (atoms.positions[index]-ads_com)/np.linalg.norm(atoms.positions[index]-ads_com)
    print(molecule_vector+atoms.positions[index])

    return molecule_vector

def rotate_mol2site_vec(atoms, atoms_idx, site_vec, mol_vec):

    # ROTATION MATRIX TAKEN FROM https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724

    import numpy as np

    theta=np.arccos(np.dot(mol_vec,site_vec))
    print(f"site and mol vectors: {site_vec}, {mol_vec}")
    axis=np.cross(mol_vec,site_vec)

    axis=axis/np.linalg.norm(axis)

    rotation_matrix=-normal_rotation_matrix(theta,axis)

    zeroed_atomco=atoms.positions-atoms.positions[atoms_idx]
    rotated_positions=np.dot(rotation_matrix, zeroed_atomco.T).T

    return rotated_positions

def normal_rotation_matrix(theta,ax):

    import numpy as np

    sinA = np.sin(theta)
    cosA = np.cos(theta)

    onemincosA = 1.0 - cosA

    rotation_matrix=np.array([[(ax[0]*ax[0]*onemincosA)+cosA ,(ax[1]*ax[0]*onemincosA)-(sinA*ax[2]) ,(ax[2]*ax[0]*onemincosA)+(sinA*ax[1])]
            ,[(ax[0]*ax[1]*onemincosA)+(sinA*ax[2]) ,(ax[1]*ax[1]*onemincosA)+cosA ,(ax[2]*ax[1]*onemincosA)-(sinA*ax[0])]
            ,[(ax[0]*ax[2]*onemincosA)-(sinA*ax[1]) ,(ax[1]*ax[2]*onemincosA)-(sinA*ax[0]) ,(ax[2]*ax[2]*onemincosA)+cosA]])

    return rotation_matrix
