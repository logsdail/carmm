def calculate_molecules(atoms,  mult=1, print_output=False):
    '''
    Returns information on the molecules in the system. Needs refinement!
    
    
    Parameters:
    atoms: Atoms object
        Input structure from which to calculate molecular information
    Mult: a multiplier for the cutoffs
        Set to 1 as default, but can be adjusted depending on application
    print_output: Boolean 
        Whether to print any information whilst working out molecules

    Returns:
    molecules: List of list of integers
        Indices of atoms involved in each molecule
        can view individual molecules via atoms[molecules[0]] ect
    '''

    from ase.neighborlist import natural_cutoffs, NeighborList
    from scipy import sparse
    cutOff = natural_cutoffs(atoms, mult=mult)
    neighborList = NeighborList(cutOff, self_interaction=False, bothways=True)
    neighborList.update(atoms)
    matrix = neighborList.get_connectivity_matrix()
    n_molecules, component_list = sparse.csgraph.connected_components(matrix)

    molecules = []
    for n in range(n_molecules):
        atomsIdxs = [i for i in range(len(component_list)) if component_list[i] == n]
        molecules.append(atomsIdxs) 
        if (print_output):
            print("The following atoms are part of molecule {}: {}".format(n, atomsIdxs))

    return molecules

def calculate_formula(atoms, mult=1):
    '''
    What it does
    Determines the chemical formula for each molecule in the system.

    Parameters:
    atoms: Atoms object
        Input structure from which to calculate molecular information.

    Returns:
    molecule_formula: List of strings.

    TO DO:
    Molecules such as CO2 are returned as C1O2. Preferable to make consistent with
    Hill's notation (CO2).
    '''

    from ase import Atoms
    from ase.formula import Formula
    molecules = calculate_molecules(atoms, mult=mult)
    molecule_formula = []
    for i in range(len(molecules)):
        # obtain chemical symbols. Ie: for CO2, they are returned as ['C', 'O', 'O'].
        symbols = [atoms[j].symbol for j in molecules[i]]
        # converts string of chemical symbols into an atoms object
        atom = Atoms(symbols)
        # uses the ASE formula functionality to convert formula into Hill notation
        chemical_formula = Formula(str(atom.symbols)).format('hill')
        # appends formula to a list containing all coordinating molecules
        molecule_formula.append(chemical_formula)

    return molecule_formula

