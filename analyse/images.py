def write_to_povray(atoms, fname):
    from ase.io import write
    write(fname, atoms)