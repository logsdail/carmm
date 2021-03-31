

    # Load Modules
    import numpy as np
    from ase.geometry.analysis import Analysis


    # loads molecules function to detect discrete molecules in atoms object
    analysis = Analysis(atoms)
    from carmm.analyse.molecules import calculate_molecules
    molecules = calculate_molecules(atoms)

    # Makes the molecules list into a *_mol atoms like object
    A = molecules[0]
    B = molecules[1]

    A_mol = atoms[A]
    B_mol = atoms[B]
    #You can view these objects separated from the rest of your system using view


    def get_lowest_distances(A_mol, B_mol,):
        '''
        Uses molecules.py to separate fn into molecules then measures the shortest distances between.
        Indented for Periodic systems

        Parameters:
         A_mol: Atoms object created by molecules
         B_mol

        Still very much a work in progress so go easy on it
                '''

    # Loops measuring the shortest distance from every atom in A to B
        measured = []
        for a in A_mol:
            bond_list_atom = []
            for b in B_mol:
                pos_diff = np.linalg.norm(a.position - b.position)
                bond_list_atom += [pos_diff]
            measured += [bond_list_atom]
        lowest_distances = [np.amin(i) for i in measured]

        return lowest_distances


print(get_lowest_distances(A_mol, 1))

view(atoms)
