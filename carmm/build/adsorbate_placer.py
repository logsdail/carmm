def place_adsorbate(atoms_site, atoms_ads, ads_idx, site_idx, bond_length, rotation=[0,0,0]):

    site_normal = find_site_normal(atoms_site, site_idx)
    adsorbate_normal = find_adsorbate_normal(atoms_mol, ads_idx)

    # Returns the zeroed coordinates about atom[ads_idx].
    atoms_site.positions = rotate_mol2site_vec(atoms_mol, ads_idx, site_normal, adsorbate_normal)

    # Rotates molecule about axis
    atoms_ads = rotate_adsorbate_site_normal(atoms_mol, atoms_site, ads_idx, site_idx, rotation)

    # Places zeroed coordinates at the bonding position
    atoms_ads.positions = atoms_site.positions + (bond_length * site_normal) + atoms_site[site_idx]

    ads_and_site = atoms_ads + atoms_site

    return ads_and_site

def rotate_adsorbate_site_normal(atoms_mol, atoms_site, ads_idx, site_idx, rotation=[0,0,0]):

    # Basic wrapper for rotation about site normal.
    # Rotations performed on a RHS axis, with rotations about site normal x, out of plane z, and orthogonal y.
    site_normal=find_site_normal(atoms_site, site_idx)
    adsorbate_normal = find_adsorbate_normal(atoms_mol, ads_idx)

    # Given the choice of two orthogonal vectors wrt the site normal is arbitrary, we choose the cross product
    # with respect to the adsorbate normal, then the cross product with this newly generated axis with the
    # site normal.
    x_axis=site_normal
    z_axis=np.cross(site_normal,adsorbate_normal)
    y_axis=np.cross(site_normal,z_axis)

    atoms_mol.rotate(rotation[0],x_axis,center=atoms_mol.positions[ads_idx])
    atoms_mol.rotate(rotation[1],y_axis,center=atoms_mol.positions[ads_idx])
    atoms_mol.rotate(rotation[2],z_axis,center=atoms_mol.positions[ads_idx])

    return atoms_mol

def find_site_normal(atoms, index):

    from carmm.analyse.neighbours import neighbours
    import numpy as np

    neighbour_atoms, shell_list = neighbours(atoms, [index], 1, verbose=False)

    vectors=atoms.positions[neighbour_atoms]-atoms.positions[index]

    site_normal=np.sum(vectors)
    site_normal=-site_normal/np.linalg.norm(site_normal)

    return site_normal

def find_adsorbate_normal(atoms, index):

    import numpy as np

    ads_com = atoms.get_center_of_mass()

    molecule_vector = (atoms[index]-ads_com)/np.linalg.norm(atoms[index]-ads_com)

    return molecule_vector

def rotate_mol2site_vec(atoms, atoms_idx, site_vec, mol_vec):

    # ROTATION MATRIX TAKEN FROM https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724

    import numpy as np

    theta=np.arccos(np.dot(site_vec,mol_vec))
    axis=np.cross(site_vec,mol_vec)

    axis=axis/np.linalg.norm(axis)

    rotation_matrix=normal_rotation_matrix(theta,axis)

    zeroed_atomco=atoms.positions-atoms.positions[atoms_idx]
    rotated_positions=np.dot(rotation_matrix, zeroed_atomco.T).T + atoms.positions[atoms_idx]

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
