'''
Short example script to measure between the center of masses of two molecules
'''



def test_analyse_center_of_masses():

    from carmm.analyse.molecules import calculate_molecules
    from carmm.analyse.planes import center_of_mass_distance
    from data.model_gen import get_example_slab as slab
    ### Traditional ASE functionality #####
    slab = slab(adsorbate=True)
    molecules = calculate_molecules(slab)
    A = molecules[0]
    B = molecules[1]
    A_mol = slab[A]
    B_mol = slab[B]
    center_of_mass_distance(A_mol, B_mol)
    #print(center_of_mass_distance(A_mol, B_mol))
    value = center_of_mass_distance(A_mol, B_mol)
    assert(1e-5 > value - 6.32082)
test_analyse_center_of_masses()