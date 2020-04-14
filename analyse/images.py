def write_to_povray(atoms, fname):
    '''
    Description

    Parameters:

    atoms: Atoms object or trajectory of individual atoms

    fname: String

    '''
    from ase.io import write
    write(fname, atoms)