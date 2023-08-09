
'''
Testing the build_mixed_oxide functionality
'''

def test_build_oxide():
    from carmm.build.mixed_oxide import mixed_oxide
    from ase.io import read, write
    import numpy as np
    orig_struc = read('/data/mixed_oxide_files/geometry.in')
    mixed_oxide_geometry = mixed_oxide(orig_struc,'Co', 'Mn',2,+2,+2, 3, 5)

    chem_sym = np.array(mixed_oxide_geometry.get_chemical_symbols())
    unique, counts = np.unique(chem_sym, return_counts=True)
    assert counts[list(unique).index('Mn')]==10

test_build_oxide()