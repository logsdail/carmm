def test_comparing_bonds_lenghts():
    from ase.io import read, write
    import numpy as np
    from ase.build import fcc111, bulk, surface, fcc110
    from carmm.analyse.bonds import comparing_bonds_lengths
    atoms1 = fcc111('Pd', (2,2,1))
    atoms2= fcc110('Pd', (2,2,1))
    atoms1.center(vacuum=20, axis=2)
    atoms1.center(vacuum=20, axis=2)
    difference = comparing_bonds_lengths(atoms1, atoms2)

    print('difference in distances is :', difference)
    print("maximium distance difference is :", difference[np.argmax(difference)])
    print("minimum distance difference is :", difference[np.argmin(difference)])
test_comparing_bonds_lenghts()