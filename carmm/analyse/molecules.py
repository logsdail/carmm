def calculate_molecules(atoms, print_output=False):
    '''
    Returns information on the molecules in the system. Needs refinement!
    
    
    Parameters:
    atoms: Atoms object
        Input structure from which to calculate molecular information 
    print_output: Boolean 
        Whether to print any information whilst working out molecules

    Returns:
    molecules: List of list of integers
        Indices of atoms involved in each molecule
    '''
    from ase.neighborlist import natural_cutoffs, NeighborList
    from scipy import sparse

    cutOff = natural_cutoffs(atoms)
    neighborList = NeighborList(cutOff, self_interaction=False, bothways=True)
    neighborList.update(atoms)
    matrix = neighborList.get_connectivity_matrix()
    n_molecules, component_list = sparse.csgraph.connected_components(matrix)

    molecules = []
    molecule_sym= []
    for n in range(n_molecules):
        atomsIdxs = [i for i in range(len(component_list)) if component_list[i] == n]
        atomsCs = [atoms[atomsIdxs[i]].symbol for i in range(len(atomsIdxs))]
        molecules.append(atomsIdxs)
        molecule_sym.append(atomsCs)
        if (print_output):
            print("Molecule {} has the following atom ids    : {}".format(n, molecules[n]))
            print("Molecule {} has the following atom symbols: {} \n".format(n, molecule_sym[n]))

    return molecules, molecule_sym
    
    #print(atoms.symbols.get_chemical_formula())
    #print(molecules[0])
    #print(type(molecules[0]))
    #print(molecules.symbols.get_chemical_formula('hill', 'empirical'))

# idx = 0
# while idx <= 251:
#    atomsIdx = component_list[idx]
#    print("There are {} molecules in the system".format(n_components))
#    print("Atom {} is part of molecule {}".format(idx, atomsIdx))
#    atomsIdxs = [i for i in range(len(component_list)) if component_list[i] == atomsIdx]
#   print("The following atoms are part of molecule {}: {}".format(atomsIdx, atomsIdxs))
#    idx = idx + 1


