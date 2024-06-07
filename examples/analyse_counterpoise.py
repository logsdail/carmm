#!/usr/bin/env python3
'''

TODO: Needs high level description

'''
def test_analyse_counterpoise():

    from os import getcwd, chdir
    from ase.io import read

    from carmm.analyse.counterpoise_onepot import counterpoise_calc
    from carmm.run.aims_calculator import get_aims_calculator
    from carmm.run.aims_path import set_aims_command

    # This is an example script for using counterpoise_calc for counterpoise (CP) correction. Please note the species
    # files in data/CO_BSSE are fake ones and default species settings are also deleted from aims.out.

    CO = read('data/CO_BSSE/C_monoxide_pbe.traj')
    examples_directory = getcwd()

    # Switch to the directory with output files, so results can be read directly without doing any actual calculation
    chdir(path='data/CO_BSSE')  # This line is only for CI-test purpose and should be deleted in actual calculation.

    # Construct the calculator
    toy_calc = get_aims_calculator(dimensions=0, xc='pbe')
    toy_calc.set(xc='pbe', spin='collinear', default_initial_moment=0.5, relativistic='atomic_zora scalar')

    # Change the species directory to current directory with fake species files
    toy_calc.set(species_dir='.')  # This line is only for CI-test purpose and should be deleted in actual calculation.

    # This function can work with lists of indices or symbols of the two parts in a binding complex for CP correction.
    # This does not work with socket calculator for now.
    # Names of species used in the CP correction can also be provided by user for clearer output.
    # Output filename would be like {a_name}_only and {a_name}_plus_ghost
    # Let's say we have A and B in this complex
    # A_only has A in the geometry of the binding complex with its own basis
    # A_plus_ghost has A in the same geometry as in the complex with B replaced by ghost atoms.

    cp_index = counterpoise_calc(CO, a_id=[1], b_id=[0], fhi_calc=toy_calc, a_name='C', b_name='O',
                                 verbose=True, dry_run=True)
    cp_symbol = counterpoise_calc(CO, a_id=['C'], b_id=['O'], fhi_calc=toy_calc, dry_run=True)

    # CP correction = A_only + B_only - A_plus_ghost - B_plus_ghost
    # This value should be added to the energy change of interest, such as adsorption energy.

    # CI-test
    assert cp_index == -8.707358006176946e-05
    assert cp_symbol == -8.707358006176946e-05

    # Check the last created geometry.in file during the calculation.
    # These three lines below are only for CI-test purpose and should be deleted in actual calculation.
    f = open("geometry.in", 'r')
    lines = f.readlines()
    assert lines[6] == "empty -0.0000000000000000 0.0000000000000000 -0.6536947973321450 C\n"

    # Return to examples directory
    chdir(path=examples_directory)  # This line is only for CI-test purpose and should be deleted in actual calculation.

test_analyse_counterpoise()

