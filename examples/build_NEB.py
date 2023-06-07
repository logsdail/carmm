#/usr/bin/env python3

'''
The functionality presented herein shows how to manipulate the ordering of atoms within Atoms objects in the context
of input preparation for nudged elastic band calculations. The initial and final image must have atoms in the correct
order to ensure that interpolation of atomic coordinates works as intended, chemical symbols must also be identical
between pairs of atoms in the initial and final Atoms object.

'''

def test_build_NEB():

    from carmm.build.neb.indices import switch_indices, switch_all_indices, sort_by_symbols
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
    updated_final = switch_indices(final, atom_to_swap, other_atom_to_swap)

    #### Assertion tests ####
    assert(not check_interpolation(initial, updated_final, 10, verbose=False, save=False))
    ########

    # To simplify the process one can also sort the indices by chemical symbols to avoid errors
    final_mismatched = switch_indices(final, 0, 1)
    try:
        assert check_interpolation(initial, final_mismatched, 10, verbose=True, save=False)
    except ValueError:
        print("Interpolation results in an error due to chemical ordering.")

    sorted_initial = sort_by_symbols(initial)
    sorted_final = sort_by_symbols(final_mismatched)

    #### Assertion tests ####
    assert (check_interpolation(sorted_initial, sorted_final, 10, verbose=True, save=False))

    # Now show how to identify indices to swap automatically, and then switch them
    #### Functionality to identify indices that would need swapping automatically
    indices_to_swap, distances = compare_structures(initial, updated_final)
    reupdated_final = switch_all_indices(updated_final, indices_to_swap)

    #### Assertion tests ####
    assert(check_interpolation(initial, reupdated_final, 10, verbose=True, save=False))



# Run the example
test_build_NEB()
