def read_dos_from_old_pickle(filename):
    '''
    Method to read in DOS data from pickle files (GPAW, outdated)

    Parameters:

    filename: String
        Name of the pickle file we are reading
    '''
    import pickle

    ef = 0
    energy = []
    dos = []

    try:
        ef, energy, dos = pickle.load(open(filename))
    except:
        print
        "A problem occurred reading the dos from ", filename
        sys.exit(1)

    return ef, energy, dos