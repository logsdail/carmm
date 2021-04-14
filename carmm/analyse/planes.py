def get_interplane_distance(atoms):
    '''
    TODO: Document (@Jack Warren)

    Returns:

    '''

    #from ase.geometry.analysis import Analysis
    # loads molecules function to detect discrete molecules in atoms object
    # AJL, Apr 2021: This unused - disabled for now. Remove?
    #analysis = Analysis(atoms)

    from carmm.analyse.molecules import calculate_molecules
    molecules = calculate_molecules(atoms)

    # Makes the molecules list into a *_mol atoms like object
    A_mol = atoms.copy()[molecules[0]]
    B_mol = atoms.copy()[molecules[1]]

    #You can view these objects separated from the rest of your system using view(A_mol)

    return get_lowest_distances(A_mol, B_mol)

def get_lowest_distances(A_mol, B_mol,):
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
