'''
Short example script to measure between the center of masses of two molecules
'''



def test_distance_between_centers_of_mass():

    from carmm.analyse.molecules import calculate_molecules
    from carmm.analyse.planes import distance_between_centers_of_mass
    from data.model_gen import get_example_slab as slab
    ### Traditional ASE functionality #####
    slab = slab(adsorbate=True)
    molecules = calculate_molecules(slab)
    A_mol = slab[molecules[0]]
    B_mol = slab[molecules[1]]
    distance_between_centers_of_mass(A_mol, B_mol)
    #print(distances_between_centers_of_mass(A_mol, B_mol))
    value = distance_between_centers_of_mass(A_mol, B_mol)
    assert(1e-5 > value - 6.32082)
test_distance_between_centers_of_mass()