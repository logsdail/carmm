def create_counterpoise_input(complex_struc=None, a_id=None, b_id=None, symbol_not_index=True):
    """
    This function creates all the geometry.in files in separate folders, assuming a binding complex
    AB.  
    Parameters:
        complex_struc: string
            This is the filename of the ase trajectory file (or any file type ase.io.read accepts)
            which stores the optimized structure of the binding complex
        a_id: list of atom symbols or atom indices for species A
        b_id: list of atom symbols or atom indices for species B
            Please use both symbols or both indices for a_id and b_id.
        symbol_not_index: bool
            This indicates whether symbols are used.
    """
    import subprocess
    from ase.io import read, write
    from ase.io.aims import write_aims
    # Let's say we have A and B in this complex
    # ?_only has A or B in the geometry of the binding complex with its own basis
    # ?_ghost has A or B in the same geometry with ghost atoms added
    dir_list = ['A_only', 'A_ghost', 'B_only', 'B_ghost']
    subprocess.run(['mkdir'] + dir_list)
    for directory in dir_list:
        binding_complex = read(complex_struc)
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
            write_aims("geometry.in", atoms=binding_complex, ghosts=ghost_list)

        subprocess.run(['mv', 'geometry.in', directory])

        # Prepare your control.in and job.sh to run FHI-aims directly before using these lines below
        # still testing
        #subprocess.run(['cp', 'control.in', directory])
        #subprocess.run(['cp', 'job.sh', directory])
        #subprocess.run(['sbatch', 'job.sh'], shell=True, cwd=directory)



