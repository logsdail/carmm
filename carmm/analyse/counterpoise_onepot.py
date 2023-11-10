def create_counterpoise_input(complex_struc=None, a_id=None, b_id=None, symbol_not_index=True):
    """
    This function creates all the geometry.in files for CP correction in separate folders,
    assuming a binding complex AB.
    Parameters:
        complex_struc: ASE Atoms
            This is the Atoms object which stores the optimized structure of the binding complex
        a_id: list
            list of atom symbols or atom indices for species A
        b_id: list
            list of atom symbols or atom indices for species B
            Please use both symbols or both indices for a_id and b_id.
        symbol_not_index: bool
            This indicates whether symbols are used.
    """
    import subprocess
    from copy import deepcopy
    from ase.io import write, aims
    # Let's say we have A and B in this complex
    # ?_only has A or B in the geometry of the binding complex with its own basis
    # ?_plus_ghost has A or B in the same geometry with ghost atoms added
    dir_list = ['A_only', 'A_plus_ghost', 'B_only', 'B_plus_ghost']
    subprocess.run(['mkdir'] + dir_list)
    for directory in dir_list:
        binding_complex = deepcopy(complex_struc)
        binding_complex.set_constraint()
        if 'A' in directory:
            ghost_id = b_id
        elif 'B' in directory:
            ghost_id = a_id

        if symbol_not_index:
            del_list = [atom.index for atom in binding_complex if atom.symbol in ghost_id]
            ghost_list = [atom.symbol in ghost_id for atom in binding_complex]
        else:
            del_list = ghost_id
            ghost_list = [atom.index in ghost_id for atom in binding_complex]

        if 'only' in directory:
            del binding_complex[del_list]
            write('geometry.in', binding_complex)
        elif 'ghost' in directory:
            aims.write_aims("geometry.in", atoms=binding_complex, ghosts=ghost_list)

        subprocess.run(['mv', 'geometry.in', directory])

    return 0



