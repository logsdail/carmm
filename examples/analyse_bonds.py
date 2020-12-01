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

test_analyse_bonds()
test_analyse_get_sorted_distances()