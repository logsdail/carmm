'''
Short example script to measure between the center of masses of two molecules
'''



def test_analyse_between_center_of_masses():
    from carmm.analyse.planes import center_of_mass_distance
    from ase.io import read
    from carmm.analyse.molecules import calculate_molecules
    ### Traditional ASE functionality #####
    output_file = "data/Dimer/aims.out"
    atoms = read(output_file)
    molecules = calculate_molecules(atoms)
    A = molecules[0]
    B = molecules[1]
    A_mol = atoms[A]
    B_mol = atoms[B]

test_analyse_between_center_of_masses()