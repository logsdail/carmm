'''
TODO: Design a process where optimal initial/final adsorbate
    arrangements are established automatically. Predict outcomes
    of symmetry operations and form a sequence before attempting
    the operations to reduce computational time.
TODO: include mirror opertion for "110" and "100" facets

'''
from carmm.build.neb.symmetry import translation, rotate_fcc, wrap_fcc
import numpy as np


def minimize_distance(initial, final, supercell_size, surface):
    '''

    Args:
        initial: Atoms object
        final: Atoms object
        supercell_size: tuple
            Tuple with 3 integers
        surface: str

    Returns:

    '''

    def evaluate_distance(atoms1, atoms2, surface):
        '''

        Args:
            atoms1:
            atoms2:
            sequence_list:

        Returns: d_start_end or None

        '''

        from scipy.spatial.distance import sqeuclidean
        from carmm.build.neb.symmetry import sort_by_xyz

        atoms1 = sort_by_xyz(atoms1, surface)

        if (atoms1.symbols == atoms2.symbols).all():
            d_start_end = sqeuclidean(atoms1.positions.flatten(),
                                      atoms2.positions.flatten())
            return d_start_end

        return None

    rotations = {"111":120, "110":180, "100":90}

    initial = wrap_fcc(initial, surface)
    final = wrap_fcc(final, surface)

    i_copy = initial.copy()

    enumeration = []

    rot_x_y_counter = np.array((0,0,0))

    d_euclidean = evaluate_distance(initial, final, surface)
    if d_euclidean:
        enumeration += [(rot_x_y_counter.copy(), d_euclidean)]

    for rot in range(int(360/rotations[surface])):
        # reset counter
        rot_x_y_counter[1], rot_x_y_counter[2] = 0, 0
        rot_x_y_counter += (1,0,0)

        initial = i_copy
        initial = rotate_fcc(initial, surface)
        i_copy = initial.copy()



        d_euclidean = evaluate_distance(initial, final, surface)
        if d_euclidean:
            enumeration += [(rot_x_y_counter.copy(), d_euclidean)]

        for x in range(supercell_size[0]):
            initial = translation(initial, axis=0, surface=surface)

            rot_x_y_counter[2] = 0
            rot_x_y_counter += (0,1,0)

            d_euclidean = evaluate_distance(initial, final, surface)
            if d_euclidean:
                enumeration += [(rot_x_y_counter.copy(), d_euclidean)]

        for y in range(supercell_size[1]):
            initial = translation(initial, axis=0, surface=surface)
            rot_x_y_counter += (0, 0, 1)

            d_euclidean = evaluate_distance(initial, final, surface)
            if d_euclidean:
                enumeration += [(rot_x_y_counter.copy(), d_euclidean)]

    enumeration = sorted(enumeration, key=lambda k: k[1])

    return enumeration


def apply_sequence(initial, final, sequence, surface):

    temp_atoms = wrap_fcc(initial, surface)
    temp_atoms2 = wrap_fcc(final, surface)

    for rot in range(sequence[0]):
        temp_atoms = rotate_fcc(temp_atoms, surface)

    for x in range(sequence[1]):
        temp_atoms = translation(temp_atoms, axis=0, surface=surface)

    for y in range(sequence[2]):
        temp_atoms = translation(temp_atoms, axis=0, surface=surface)

    return temp_atoms, temp_atoms2

