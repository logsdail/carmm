# Authors: Owain Beynon, Igor Kowalec
def neighbours(atoms, centre, shell, cutoff=None, verbose=False):
    ''' Returns a list of indices of atomic neighbors from a central atom

    Parameters:
    atoms : Atoms object
        Input structure to count neighbours
    centre : list of integers
        Indices of atom(s) to start counting from
    shell : Integer
        Size of the nearest neighbour shell, 1st neighbours, 2nd neighbours ....
    cutoff: list of floats or None
        Bond length cutoff distance in Angstrom can be set for each atom individually
        The list must contain exactly len(atoms) floats. If None, natural_cutoffs
        are used by default.
    verbose: boolean
        If True, information about the cutoffs and selected neighbors is printed

    Returns:
        List of all atoms within specified neighbour distances
        List of lists containing indices of atoms that make 0th, 1st (...) shell
    '''

    from ase.neighborlist import natural_cutoffs, NeighborList

    if not cutoff:
        # Choose a cutoff to determine max bond distance, otherwise use natural cutoff
        cutoff = natural_cutoffs(atoms)
        if verbose:
            print("Default bond cutoffs selected:", set([(atoms[i].symbol, cutoff[i]) for i in range(len(atoms))]))

    # Set object storing all neighbour information
    all_neighbours = set(centre)
    # List of lists storing each neighbour shell requested
    shell_list = [centre]

    # Creates an empty list an appends atom indices whose distances are
    # x amount nearest neighbours away from centre
    for this_shell in range(shell):
        # keep new neighbor indices in a set to avoid duplicates
        new_neighbors = set()
        for index in all_neighbours:
            # find neighbors based on cutoff and connectivity matrix
            nl = NeighborList(cutoff, self_interaction=False, bothways=True)
            nl.update(atoms)
            indices = nl.get_neighbors(index)[0]
            for i in indices:
                new_neighbors.add(i)

        shell_list += [[i for i in new_neighbors if i not in all_neighbours]]
        for i in new_neighbors:
            all_neighbours.add(i)

        if verbose:
            for shell in range(len(shell_list)):
                print("Shell", shell, "contains atoms with indices:", shell_list[shell])

        all_neighbours = list(all_neighbours)
        atoms_copy = atoms.copy()
        selection = atoms_copy[all_neighbours]

    return all_neighbours, shell_list, selection


# Authors: Igor Kowalec, Lara Kabalan, Jack Warren
#
def surface_coordination(atoms, cutoff=None, verbose=True):
    '''
    This function allows to extract the following data from the
    supplied Atoms object in the form of a dictionary:
    Per atom:
    - number of neighbouring atoms of specific chemical symbol

    Per atomic layer (based on tag):
    - average coordination number for all M-M combinations of all chemical species,
        e.g. for CuAu alloy slab an average number of Cu-Cu, Cu-Au, Au-Cu and Au-Au
            bonds with surrounding atoms i.e. Cu-Cu_neighbors_per_layer etc.
    - concentration of atoms per layer


    Parameters:
        atoms: Atoms object
            Surface slab containing tagged atomic layers
            e.g. using carmm.build.neb.symmetry.sort_z
        cutoff: list of floats or None
             Bond length cutoff distance in Angstrom can be set for each atom individually
             The list must contain exactly len(atoms) floats. If None, natural_cutoffs
             are used by default.
        verbose: boolean
            If True, analysed data that is contained in the dictionary will be printed
            as a table
            TODO: make table neater
    Returns:
        dict_CN, dict_surf_CN
            TODO: proper description of dictionary structure, writing csv files'''
    from ase.neighborlist import natural_cutoffs, NeighborList
    import numpy as np
    from collections import Counter
    from itertools import product

    dict_CN = {}
    if not cutoff:
        # Choose a cutoff to determine max bond distance
        cutoff = natural_cutoffs(atoms)
        if verbose:
            print("Default bond cutoffs selected:", set([(atoms[i].symbol, cutoff[i]) for i in range(len(atoms))]))

    # Create a symbols set based on all species present in the atoms
    symbols_set = sorted(list(set(atoms.symbols)))

    for l in set(atoms.get_tags()):
        index= [i.index for i in atoms if i.tag == l]
        for i in index:
            nl = NeighborList(cutoff, self_interaction=False, bothways=True)
            nl.update(atoms)
            indices, offsets = nl.get_neighbors(i)
            # avoid duplicate atomic indices
            indices_no_self = np.array(list(set([j for j in indices])))

            cn_numbers = Counter(atoms[indices_no_self].symbols)

            # create a dictionary of relevant values for cn calculations
            dict_CN[i] = {"symbol":atoms[i].symbol, 'index':i, "layer":l}
            for k in symbols_set:
                dict_CN[i].update({k+"_neighbors":cn_numbers[k]})

    # check all combinations of atomic symbols
    pr = product(symbols_set, repeat=2)

    dict_surf_CN= {}
    for layer in set(atoms.get_tags()):
        dict_surf_CN[layer] = {}
        dict_surf_CN[layer].update({"layer":layer})

        # extract data for all combinations of neighbors and put at the end of the dictionary
        for p in product(symbols_set, repeat=2):
            M_M_cn = [dict_CN[x][p[0] + "_neighbors"] for x in dict_CN if\
                 dict_CN[x]["layer"] == layer and dict_CN[x]['symbol'] == p[1]]
            if not M_M_cn == []:
                M_M_avg_cn = np.average(M_M_cn)
            else:
                M_M_avg_cn = "0" # as the concentration of one of the metals is zero in the layer

            dict_surf_CN[layer].update({p[0] + "_neighboring_w_" + p[1]: M_M_avg_cn})

        # Check concentrations of atoms per layer
        for symbol in symbols_set:
            atoms_per_layer = len([x for x in dict_CN if dict_CN[x]["layer"] == layer])
            if atoms_per_layer > 0:
                symbol_count_per_layer = len([symbol for x in dict_CN if dict_CN[x]["layer"] == layer and dict_CN[x]['symbol'] == symbol])
                symbol_concentration_per_layer = symbol_count_per_layer / atoms_per_layer
            else:
                symbol_concentration_per_layer = 0
            dict_surf_CN[layer].update({symbol+"_concentration_per_layer":symbol_concentration_per_layer})

    # Preparation for verbose
    cn_layer_list_dict = [dict_surf_CN[i] for i in range(len(dict_surf_CN))]
    cn_layer_dict_keys = [dict_surf_CN[i].keys() for i in range(len(dict_surf_CN))][0]

    # Optional text output
    if verbose:
        print([key for key in cn_layer_dict_keys])
        for i in range(len(dict_surf_CN)):
            print([cn_layer_list_dict[i][key] for key in cn_layer_dict_keys])


    return dict_CN, dict_surf_CN
