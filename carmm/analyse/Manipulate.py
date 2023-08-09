def Crystal_rotations(atoms, ):
    import numpy as np
    from ase.build import rotate
    from ase.io import read
    from ase.visualize import view
    from carmm.analyse.molecules import calculate_molecules
    from carmm.analyse.planes import establish_planes
    from sklearn.decomposition import PCA

    atoms = read("filename")
    molecules = calculate_molecules(atoms)
    view(atoms)
    A_mol = atoms[molecules[0]]
    B_mol = atoms[molecules[1]]
    C_mol = atoms[molecules[2]]
    D_mol = atoms[molecules[3]]
    E_mol = atoms[molecules[4]]
    #print(molecules)
    #view(A_mol)
    #view(B_mol)
    #view(C_mol)
    #view(D_mol)
    #view(E_mol)
    m1 = A_mol.positions
    #m2 = B_mol.get_positions(B_mol)
    pca = PCA(n_components=2)
    pca.fit(m1)

    #print (pca.components_)

    Atom1 = m1[102]
    Atom2 = m1[103]
    Atom3 = m1[104]

    Vector1 = Atom1 - Atom2
    Vector2 = Atom2 - Atom3

    cross = np.cross(*pca.components_)
    crossBN = np.cross(Vector2,Vector1)
    #Dot1 = -np.sum(Atom1*Vector) # dot product
    #print(Vector1,Vector2,cross)
    #z = (-Vector[0] - Vector[1]- Dot)*1./Vector[2]
    CM_A = A_mol.get_center_of_mass()
    CM_C = D_mol.get_center_of_mass()
    cross *= np.dot(CM_C - CM_A, cross)
    crossBN *=np.dot(CM_C - CM_A, crossBN)
    C_mol.translate(-CM_C)
    #print(Vector)
    C_mol.rotate(30,crossBN,'COP',rotate_cell=False)
    C_mol.translate(CM_C)
    recom = A_mol + C_mol
    view(recom)
    Fin = recom + B_mol + D_mol + E_mol
    view(Fin)

    #print(CM_B,CM_A)