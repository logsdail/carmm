import numpy as np


def get_interplane_distance(atoms):
    '''
    TODO: Document (@Jack Warren)

    Returns:

    '''

    # from ase.geometry.analysis import Analysis
    # loads molecules function to detect discrete molecules in atoms object
    # AJL, Apr 2021: This unused - disabled for now. Remove?
    # analysis = Analysis(atoms)

    from carmm.analyse.molecules import calculate_molecules
    molecules = calculate_molecules(atoms)

    # Makes the molecules list into a *_mol atoms like object
    A_mol = atoms.copy()[molecules[0]]
    B_mol = atoms.copy()[molecules[1]]

    # You can view these objects separated from the rest of your system using view(A_mol)

    return get_lowest_distances(A_mol, B_mol)


def get_lowest_distances(A_mol, B_mol, ):
    '''
    Uses molecules.py to separate fn into molecules then measures the shortest distances between.
    Indented for Periodic systems

    Parameters:
    A_mol: Atoms object created by molecules
    B_mol: TODO: Complete

    Returns:

    Still very much a work in progress so go easy on it
    '''
    import numpy as np

    # Loops measuring the shortest distance from every atom in A to B
    measured = []
    for a in A_mol:
        bond_list_atom = []
        for b in B_mol:
            pos_diff = np.linalg.norm(a.position - b.position)
            bond_list_atom += [pos_diff]
        measured += [bond_list_atom]
    lowest_distances = [np.amin(i) for i in measured]

    # TODO: Document the return - this is a list of shortest distances from each atom in A
    # TODO: What about information on the atoms that constitute the lowest distance?
    return lowest_distances


def center_of_mass_distance(A_mol, B_mol, ):
    '''

    :param A_mol:
    :param B_mol:
    :return:  distance between 2 centers of mass
    '''
    # Gets the Center of mass for the molecules and simply measures between
    CM_A = A_mol.get_center_of_mass()
    CM_B = B_mol.get_center_of_mass()
    pos_diff = np.linalg.norm(CM_A - CM_B)
    return pos_diff
    print(pos_diff)


def establish_planes(Atom1, Atom2, Atom3,):
    '''
    Using 3 points to calculate 3 vectors and establish a plane...

    :param Atom1: (x,y,z coordinates of Atom1)
    :param Atom2:
    :param Atom3:
    :return: Graphical Representation of planes
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    Atom1 = np.array([5,5,5])
    Atom2 = np.array([0,-4,0])
    Atom3 = np.array([0,0,1])

    Vector1 = Atom1-Atom2
    Vector2 = Atom2-Atom3
    Vector3 = Atom1-Atom3

# a plane is a*x+b*y+c*z+d=0
# [a,b,c] is the normal. Thus, we have to calculate d and we're set

    d1 = -np.sum(Atom1*Vector1) # dot product
    d2 = -np.sum(Atom2*Vector2) # dot product
    d3 = -np.sum(Atom3*Vector3) # dot product

# create x,y
    xx, yy = np.meshgrid(range(30), range(30))

# calculate corresponding z
    z1 = (-Vector1[0]*xx - Vector1[1]*yy - d1)*1./Vector1[2]
    z2 = (-Vector2[0]*xx - Vector2[1]*yy - d2)*1./Vector2[2]
    z3 = (-Vector3[0]*xx - Vector3[1]*yy - d3)*1./Vector3[2]

# plot the surface
    plt3d = plt.figure().gca(projection='3d')
    plt3d.plot_surface(z1,z2,z3, color='blue')
#    plt3d.plot_surface(xx,yy,z2, color='red')
#    plt3d.plot_surface(xx,yy,z3, color='green')
    plt.show()