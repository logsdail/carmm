def find_site_normal(atoms, index):

    from carmm.analyse.neighbours import neighbours
    import numpy as np

    neighbour_atoms, shell_list = neighbours(atoms, [index], 1, verbose=False)

    vectors=atoms.positions[neighbour_atoms]-atoms.positions[index]

    site_normal=np.sum(vectors)
    site_normal=-site_normal/np.linalg.norm(site_normal)

    return site_normal

def molecule_normal(atoms, index):

    import numpy as np

    ads_com = atoms.get_center_of_mass()

    molecule_vector = (atoms[index]-ads_com)/np.linalg.norm(atoms[index]-ads_com)

    return molecule_vector

def rotate_mol2site_vec(site_vec, mol_vec):

    import numpy as np

    np.dot(site_vec,mol_vec)
    np.cross(site_vec,mol_vec)



def place_adsorbate_site_vec(atoms, vec):


define a center of rotation

define sanity test for close species

EXTRA: define complex geometries to snap to (eg. hexagons)

EXTRA: adjustable plot with sliders

define species geometry by neighbours

import numpy as np

#ROTATION MATRIX TAKEN FROM https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724

A=np.array([0.707,0.707,0.0])
B=np.array([0.0,1.0,0.0])

ax = np.cross(A,B)

ax = ax/np.linalg.norm(ax)
dot = np.dot(A,B)

angle=np.arccos(dot)

rotatedA=np.dot(rotation.T,B)

def normal_rotation_matrix(theta):

    sinA = np.sin(theta)
    cosA = np.cos(theta)

    onemincosA = 1.0 - cosA

    rotation_matrix=np.array([[(ax[0]*ax[0]*onemincosA)+cosA ,(ax[1]*ax[0]*onemincosA)-(sinA*ax[2]) ,(ax[2]*ax[0]*onemincosA)+(sinA*ax[1])]
            ,[(ax[0]*ax[1]*onemincosA)+(sinA*ax[2]) ,(ax[1]*ax[1]*onemincosA)+cosA ,(ax[2]*ax[1]*onemincosA)-(sinA*ax[0])]
            ,[(ax[0]*ax[2]*onemincosA)-(sinA*ax[1]) ,(ax[1]*ax[2]*onemincosA)-(sinA*ax[0]) ,(ax[2]*ax[2]*onemincosA)+cosA]])

    return rotation_matrix
