#!/usr/bin/env python3

'''
This gives example usage when looking at how to get bond information from an Atoms object

This comes in useful when analysing input/output structures
'''

def test_analyse_bonds():

    from carmm.analyse.bonds import analyse_all_bonds

    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)
    # Deform the z-coordinate for C so the Au-C bond length is too short, to test deformation
    slab.positions[-3][2] -= 1.5
    #########

    abnormal_count, abnormal_bond_names = analyse_all_bonds(slab, abnormal=True)

    assert(len(abnormal_count) == 1)

def test_analyse_get_sorted_distances():
    from carmm.analyse.bonds import get_sorted_distances

    #Build a model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    distances = get_sorted_distances(slab)
    # Check values are stil lthe same and ordering is correct
    assert(len(distances) == 210)
    assert(1e-5 > abs(distances[0] - 1.178657))
    assert(1e-5 > abs(distances[-1] - 7.132005))

    # Test plot
    from carmm.analyse.distribution_functions import plot_distribution_function
    plt = plot_distribution_function(distances, title='Radial Distribution Function')
    #plt.show()

def test_analyse_comparing_bond_lengths():
    import numpy as np
    from ase.build import fcc111, fcc110
    from carmm.analyse.bonds import comparing_bonds_lengths

    # Build some models
    atoms1 = fcc111('Pd', (2, 2, 1), vacuum=20.0)
    atoms2 = fcc110('Pd', (2, 2, 1), vacuum=20.0)

    difference = comparing_bonds_lengths(atoms1, atoms2)

    # Assertion test
    assert(1e-5 > difference[np.argmax(difference)] - 2.013612)

    #print('difference in distances is :', difference)
    #print("maximium distance difference is :", difference[np.argmax(difference)])
    #print("minimum distance difference is :", difference[np.argmin(difference)])

def test_analyse_chelation():
    ## Initialises modules
    from carmm.analyse.bonds import analyse_chelation
    from ase import Atoms
    import numpy as np
    ## Builds a test hexahydrate complex to investigate coordination type.
    p = np.array(
        [[0, 0, 0],
        [0, 0,  2.0], [0,  0.7,  2.7], [0, -0.7,  2.7],
        [0, 0, -2.0], [0,  0.7, -2.7], [0, -0.7, -2.7],
        [0,  2.0, 0], [0,  2.7,  0.7], [0,  2.7, -0.7],
        [0, -2.0, 0], [0, -2.7, -0.7], [0, -2.7,  0.7],
        [ 2.0, 0, 0], [ 2.7, -0.7, 0], [ 2.7,  0.7, 0],
        [-2.0, 0, 0], [-2.7,  0.7, 0], [-2.7, -0.7, 0]])
    atoms = Atoms('Mn(OH2)6', positions=p)
    ## calculates chelation types
    ligands = analyse_chelation(atoms=atoms, metal='Mn', ligand_atoms=['O'], mult=1.5)
    ## assertion check to verify results
    assert(ligands.get("complex") == "Mn(κ1-H2O)6")

test_analyse_bonds()
test_analyse_get_sorted_distances()
test_analyse_comparing_bond_lengths()
test_analyse_chelation()
