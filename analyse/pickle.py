def read_dos_from_files(files):
    ef, energy, dos = read_dos_from_file(filename)

    print
    "Ef: ", ef
    for i in range(len(energy)):
        print
        energy[i], ",", dos[i]


def read_dos_from_old_pickle(filename):
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