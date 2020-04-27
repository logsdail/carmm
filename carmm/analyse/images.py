def write_to_povray(atoms, fname):
    '''
    Saves the incoming atoms object to a povray file

    Parameters:

    atoms: Atoms object or trajectory of individual atoms
        Structure to write as Povray
    fname: String
        Output file name, no file extension necessary.
    '''
    from ase.io import write
    write(fname+'.pov', atoms, run_povray=True)