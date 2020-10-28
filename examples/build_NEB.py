#/usr/bin/env python3

'''
TODO: Description Needed

'''

def test_build_NEB():

    from carmm.build.neb.indices import switch_indices, switch_all_indices
    from carmm.build.neb.interpolation import check_interpolation
    from carmm.analyse.bonds import compare_structures
    from data.model_gen import get_example_adsorbate as co2

    initial = co2()
    final = co2()

    # Check whether the initial pathway is valid
    #### Assertion tests ####
    assert(check_interpolation(initial, final, 10, verbose=False, save=False))
    ########

    # Swap two atoms to show how this is done
    atom_to_swap = 1
    other_atom_to_swap = 2
    #updated_final = switch_indices('final.traj', atom_to_swap, other_atom_to_swap):
    updated_final = switch_indices(final, atom_to_swap, other_atom_to_swap)

    #### Assertion tests ####
    assert(not check_interpolation(initial, updated_final, 10, verbose=False, save=False))
    ########

    # Now show how to identify indices to swap automatically, and then switch them
    #### Functionality to identify indices that would need swapping automatically
    indices_to_swap, distances = compare_structures(initial, updated_final)
    reupdated_final = switch_all_indices(updated_final, indices_to_swap)

    #### Assertion tests ####
    assert(check_interpolation(initial, reupdated_final, 10, verbose=True, save=False))

# Run the example
test_build_NEB()
