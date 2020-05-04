'''This file is work in progress'''


def translation(model, axis=0, surface="111"):
    '''
    Performs a translaton of the model by manipulation of the unit cell,
    maintaining the optimised geometry and forces. After translation original
    Atoms object is permanently changed.

    TODO:
    - FOR NOW requires surface to have tags for layers of atoms in Z-direction
    - functionality beyond FCC? or higher index?

    Parameters:
    # TODO: This routine changes model in place (inadvertently).
            Should it work with a copy of model for the "new model", to prevent editing of the old model?
    model: Atoms object
        XXX
    a: float
        lattice parameter used in the model
    axis: integer
        Choice - 0, 1 representing axis x, y
    surface: string
        FCC surface - so far supports "111", "110", "100"
    '''

    import numpy as np
    import math

    # Retrieve constraints, calculator from the model for later
    constraint = model._get_constraints()
    prev_calc = model.get_calculator()

    '''Section on variables'''
    indices_to_move = []
    # Extraction of lattice parameter from bottom layers
    shift_dist = get_lattice_constant(model)

    if surface == "110":
        shift_x = shift_dist
        shift_y = shift_dist
    elif surface == "111":
        shift_x = shift_dist * math.cos(math.radians(60))
        shift_y = shift_dist * math.cos(math.radians(30))
    elif surface == "100":
        shift_x = shift_dist * math.cos(math.radians(30))
        shift_y = shift_dist

    '''Section on moving atoms'''
    if axis == 0:
        # TODO: get rid of the rotation, it is there because math is easier
        #       would require editing sort_y_tag function
        # align atoms perpendicular to x-axis
        if surface == "111":
            model.rotate(30, 'z', rotate_cell=True)

        positions = model.get_positions()
        count_iter = 0
        for coordinate in positions:
            count_iter = count_iter + 1
            # Include tolerance
            if coordinate[0] < shift_x*1.02:
                indices_to_move = indices_to_move + [count_iter - 1]

        if surface == "111":
            model.rotate(-30, 'z', rotate_cell=True)

        # rpt cell, take from one side and add temp, exchange change positions
        temp_model = model.repeat((2, 1, 1))
        temp_model.set_constraint()
        model.set_constraint()
        temp_indices = np.array(indices_to_move) + len(model.get_tags())
        # Merge list of indices to switch
        index_pairs = zip(indices_to_move, temp_indices)

        for pairs in index_pairs:
            model[pairs[0]].position = temp_model[pairs[1]].position

        del temp_model

        if surface == "111":
            model.set_positions(
                            np.array(model.get_positions()) - np.array(
                                (shift_dist, 0, 0)))
        elif surface == "100":
            model.set_positions(
                            np.array(model.get_positions()) - np.array(
                                (shift_dist, 0, 0)))
        elif surface == "110":
            model.set_positions(
                            np.array(model.get_positions()) - np.array(
                                (shift_dist*2/(2**1/2), 0, 0)))
    elif axis == 1:
        # identify atoms to be moved
        positions = model.get_positions()
        count_iter = 0
        for coordinate in positions:
            count_iter = count_iter + 1
            # Include some tolerance for surface atoms
            if coordinate[1] < shift_y*0.98:
                indices_to_move = indices_to_move + [count_iter - 1]

        # rpt cell, take from one side and add temp, exchange change positions
        temp_model = model.repeat((1, 2, 1))
        temp_model.set_constraint()
        model.set_constraint()
        temp_indices = np.array(indices_to_move) + len(model.get_tags())
        # Merge list of indices to switch
        index_pairs = zip(indices_to_move, temp_indices)

        for pairs in index_pairs:
            model[pairs[0]].position = temp_model[pairs[1]].position

        del temp_model

        if surface == "111":
            model.set_positions(
                            np.array(model.get_positions()) - np.array(
                                (shift_x, shift_y, 0)))
        elif surface == "100":
            model.set_positions(
                            np.array(model.get_positions()) - np.array(
                                (0, shift_y, 0)))
        elif surface == "110":
            model.set_positions(
                            np.array(model.get_positions()) - np.array(
                                (0, shift_y, 0)))

    else:
        raise ValueError
        print("Axis index must be an integer - 0 or 1.")


    '''Section on index correction'''
    model = sort_by_xyz(model, surface)

    # TODO: Zero should be based on similarity to the surface atom. Does not
    #   work as intended for odd number of layers, might not be required with
    #   get_a in place
    # No need for these imports as the functions are in this file!
    # from carmm.build.neb import get_zero_from_constrained_atoms
    zero_index = get_zero_from_constrained_atoms(model)
    import copy
    zero_x = copy.deepcopy((model[zero_index].position[0]))
    zero_y = copy.deepcopy((model[zero_index].position[1]))
    # Reset to zero
    for i in reversed([atom.index for atom in model]):
        model.positions[i][0] = (model.positions[i][0] - zero_x)
        model.positions[i][1] = (model.positions[i][1] - zero_y)

    # Retain calculator information and constraint
    if model.get_calculator() is not None:
        prev_calc.atoms = model
    model.set_calculator(prev_calc)
    model.set_constraint(constraint)

    model = sort_by_xyz(model, surface)

    return model


def get_lattice_constant(model):
    '''
    Retrieve minimum M-M distance from bottom layer of a slab

    Parameters:
    model: Atoms Object
        FCC surface model that has atom tags associated with layers
    '''
    import numpy as np
    max_tag = np.amax(list(set(model.get_tags())))

    # need to convert array into list
    bottom_layer = ([atom.index for atom in model if atom.tag == max_tag])

    dist_matrix = model[bottom_layer].get_all_distances()
    a = np.amin(np.setdiff1d(dist_matrix, np.array(0.0)))
    return a


def sort_by_xyz(model, surface):
    '''
    Sorting indices by xyz coordinates for periodic surface models.
    Returns Atoms object with atom indices sorted.

    Parameters:
    model: Atoms object
        periodic surface model, so far FCC100, 111, 110 supported
    surface: string
        Face centered cubic low index surfce - "111", "110" or "110"
    '''
    import numpy as np

    if model.get_calculator() is not None:
        # Retrieve forces for forces array adjustments
        calc = model.get_calculator()
        calc_results = calc.results
        if "forces" in calc_results:
            f = calc_results["forces"]
        else:
            f = []

    # Sorting mechanism for Y tags
    def sort_y_tag(xyz, surface):
        import numpy as np
        # TODO: shortest_AB could be extracted from bottom fixed layer for
        #    consistency, but poses problems if layers are not the same
        if surface == "111":
            y_tag = np.int(np.round((2*xyz[1])/(
                        np.sin(np.radians(60))*np.amin(shortest_ABs))))
        if surface == "110":
            y_tag = np.int(np.round((2*xyz[1])/(np.amin(shortest_ABs))))

        if surface == "100":
            y_tag = np.int(np.round((2*xyz[1])/(np.amin(shortest_ABs))))
        return y_tag

    # sort z direction by tags
    # TODO: identify z-tags like y-tags above, for models not created in ASE
    # reverse, need to start from bottom, ie. last tag
    tags = list(set(model.get_tags()))
    tags = sorted(tags, reverse=True)
    index_by_xyz = []

    for tag in tags:
        index_by_tag = [atom.index for atom in model if atom.tag == tag]
        xyz = model[index_by_tag].get_positions()

        if not tag == 0:
            all_AB = model[index_by_tag].get_all_distances()
            shortest_ABs = []
            # Figure minimum bond lengths in this set of atoms
            for distance in all_AB:
                # Need to remove 0.0 from array before finding min value
                distance = np.amin(np.setdiff1d(distance, np.array(0.0)))
                shortest_ABs += [distance]

            sorted_xyz = sorted(
                xyz, key=lambda k: [sort_y_tag(k, surface), k[0]])  # Y-tag
        else:
            sorted_xyz = xyz

        # Identify sequence that will correctly reorder Atoms based on xyz
        for coordinate in sorted_xyz:
            for index in index_by_tag:
                if (coordinate == model[index].position).all():
                    index_by_xyz += [index]

    model = model[index_by_xyz]
    # Ensure function works if force information empty and rearrange
    # TODO: CHECK IF THIS IS NECESSARY
    if model.get_calculator() is not None:
        f = f[index_by_xyz]

    return model


def get_zero_from_constrained_atoms(model, mirror=False):
    '''
    Function retrieving and index of an atom in the bottom layer of surface
    supplied. Based on the assumption that model:
    i) contains constrained bottom layers
    ii) contains tags associated with layers in Z-direction
    If no contraints - align with the top layer atom.

    Parameters:
    model: Atoms object
        FCC surface model
    mirror: boolean
        If True chooses a different layer to account for mirroring of the model

    TODO: What if no tags associated with layers in z-direction?
    '''
    import numpy as np
    zero_index = None
    tag = 1
    if mirror is True:
        tag = 2
    first_layer_atom = [atom.index for atom in model if atom.tag == tag][0]
    if not model._get_constraints() == []:
        prev_const = model._get_constraints()
        # retrieve indices of the fixed atoms
        constrained_indices = np.ndarray.tolist(prev_const[0].index)
        # allow room for surface layer distortion when matching
        eps = 0.3

        sort_y = [atom.index for atom in model[constrained_indices] if
            eps > atom.position[1] - model[first_layer_atom].position[1] > -eps]

        zero_index = [atom.index for atom in model[sort_y] if
            eps > atom.position[0] - model[first_layer_atom].position[0] > -eps][0]

        zero_index = sort_y[zero_index]

    # if no constrained layers detected - return index of 1st atom in top layer
    if zero_index is None:
        zero_index = first_layer_atom
    else:
        pass

    return zero_index


def check_for_negative_positions(model, axis, surf=None):
    '''
    Check if negative positions persist after a symmetry operation
    and if yes - perform a translation until all atoms are wrapped into
    the unit cell

    Parameters:
    model: Atoms object
    axis: 0 or 1 # for "x" and "y"
    surf: str
        "111", "110", "100"
        fcc surface "111" has a special operation for axis = 0

    '''

    still_negative = False

    # align  atoms in y for x-coordinate check
    if surf == "111":
        if axis == 0:
            model.rotate(30, 'z', rotate_cell=True)
        else:
            pass
    else:
        pass

    # check if there are negative coordinates
    for i in [atom.index for atom in model]:
        if model[i].position[axis] < -1.0:
            still_negative = True

    if surf == "111":
        if axis == 0:
            model.rotate(-30, 'z', rotate_cell=True)
        else:
            pass
    else:
        pass

    x = 0
    while still_negative is True:
        x += 1
        model = translation(model, axis, surf)
        check = 0

        if surf == "111":
            if axis == 0:
                model.rotate(30, 'z', rotate_cell=True)
            else:
                pass
        else:
            pass

        # check for negative coordinates outside of unit cell in axis
        for i in [atom.index for atom in model]:
            if model[i].position[axis] < -1.0:
                check += 1

        if surf == "111":
            if axis == 0:
                model.rotate(-30, 'z', rotate_cell=True)
            else:
                pass
        else:
            pass

        # if no negative coordinates - stop
        if check == 0:
            still_negative = False
        else:
            pass

    return model


def mirror(model, center_index, plane="y", surf="111"):
    '''
    Function that returns a mirror image of an FCC model in the x or y axis
    with respect to an adsorbate and atom shifts the surface atoms accordingly.

    Parameters:
    model: Atoms object
        # TODO: This routine I think inadvertently edits the incoming model
        #       Should we make a copy as the model comes in, and work exclusively with that?
        periodic FCC surface model with an adsorbate. Tags required for surface
        layers.
    plane: string
        "x" or "y"
    center_index: integer
        Index of an atom with respect to which the image will be
        mirrored.
    surf: string
        Surface type. Supported "111", "110", "100" of the FCC lattice
    '''
    if plane == "x":
        axis = 0
    elif plane == "y":
        axis = 1
    else:
        print("Wrong axis chosen. Please choose 'x' or 'y'")

    # Retrieve calculator information
    if model.get_calculator() is not None:
        prev_calc = model.get_calculator()

    # Save initial center atom position
    # Hypothesis: center_atom_position is a shallow reference that does not
    #    survive other operations unchanged.
    # Solution: a deepcopy required to retain functionality? Confirmed.
    import copy
    center_atom_position = copy.deepcopy(model[center_index].position)

    if axis == 0:
        translate = model.positions[center_index][1]
        for i in [atom.index for atom in model]:
            model.positions[i][1] = (-model.positions[i][1] + (2*translate))
    if axis == 1:
        translate = model.positions[center_index][0]
        for i in [atom.index for atom in model]:
            model.positions[i][0] = (-model.positions[i][0] + (2*translate))

    # Sort after changing the model prior to alignment
    # from carmm.build.neb import sort_by_xyz
    model = sort_by_xyz(model, surf)

    # align to position zero - constrained atoms aligned with surface atoms
    if plane == "y":
        zero_index = get_zero_from_constrained_atoms(model, mirror=True)
    elif plane == "x":
        zero_index = get_zero_from_constrained_atoms(model)

    zero_x = (model[zero_index].position[0])
    zero_y = (model[zero_index].position[1])

    for i in reversed([atom.index for atom in model]):
        model.positions[i][0] = (model.positions[i][0] - zero_x - model.get_cell_lengths_and_angles()[0])
        model.positions[i][1] = (model.positions[i][1] - zero_y - model.get_cell_lengths_and_angles()[1])

    model = sort_by_xyz(model, surf)

    # Return to around the position of the center_index
    # Force translations to remove inconsistencies

    model = check_for_negative_positions(model, 1, surf)
    model = check_for_negative_positions(model, 0, surf)

    a = get_lattice_constant(model)
    current_pos = model[center_index].position

    # Safety break
    x = 0
    # Align in y-direction
    a = get_lattice_constant(model)
    while not (a/2 > (current_pos[1] - center_atom_position[1]) > -a/2):
        x = x+1
        current_pos = model[center_index].position
        model = translation(model, axis=1, surface=surf)
        if x > 9:
            break

    # Align in x-direction
    x = 0
    while not a/2 > (current_pos[0] - center_atom_position[0]) > -a/2:
        x = x+1
        current_pos = model[center_index].position
        model = translation(model, axis=0, surface=surf)
        if x > 9:
            break

    model = sort_by_xyz(model, surf)

    if model.get_calculator() is not None:
        prev_calc.atoms = model
        model.set_calculator(prev_calc)

    return model


def rotate_fcc(model, center_index, surf):
    '''
    Rotate FCC cell with respect to an atom by allowed increments, ie.
    "111" - 120 degrees
    "100" - 90 degrees
    "110" - 180 degrees

    Parameters:
    model:in Atoms object
        # TODO: This inadvertently edits the old model.
        #       Should we create a copy at the very start to prevent editing?
        FCC low index surface, (111), (110) or (100)
    center_index: int
        index of an atom as the center of the rotation operation
    surf: str
        "111", "110", "100" allowed
    '''
    #from carmm.build.neb import sort_by_xyz
    import numpy as np
    import copy

    if surf == "111":
        degrees = -120
    elif surf == "100":
        degrees = -90
    elif surf == "110":
        degrees = -180
    else:
        pass
        # TODO Error message

    # align to position zero
    # WIP
    if not model._get_constraints() == []:
        prev_const = model._get_constraints()
        # retrieve indices of the fixed atoms
        constrained_indices = np.ndarray.tolist(prev_const[0].index)
        eps = 1e-8

        index_const_at_zero = None
    else:
       prev_const = None
       index_const_at_zero = None

    # Remove constraints for rotation operation
    model.set_constraint()
    # this variable does not survive other operations unchanged, hence deepcopy
    center_atom_position = copy.deepcopy(model[center_index].position)

    model.rotate(degrees, 'z', center=center_atom_position)
    model = sort_by_xyz(model, surf)

    if index_const_at_zero is not None:
        zero_x = index_const_at_zero.position[0]
        zero_y = index_const_at_zero.position[1]
    else:
        zero_x = (model[[atom.index for atom in model if atom.tag == 1][0]].position[0])
        zero_y = (model[[atom.index for atom in model if atom.tag == 1][0]].position[1])
    # Reset to zero
    for i in reversed([atom.index for atom in model]):
        model.positions[i][0] = (model.positions[i][0] - zero_x)
        model.positions[i][1] = (model.positions[i][1] - zero_y)

    # Wrap function works better after alignment for fcc111 - could be
    # coincidental, but seems to work for the example provided
    model.wrap()

    # Return to around the position of the center_index
    # Force translations to remove inconsistencies
    a = get_lattice_constant(model)
    current_pos = model[center_index].position

    if a/2 > (current_pos[1] - center_atom_position[1]) > -a/2:
        model = translation(model, axis=1, surface=surf)
        current_pos = model[center_index].position
    if a/2 > (current_pos[0] - center_atom_position[0]) > -a/2:
        model = translation(model, axis=0, surface=surf)
        current_pos = model[center_index].position

    # Safety break
    x = 0
        # Align in y-directio
    while not (a/2 > (current_pos[1] - center_atom_position[1]) > -a/2):
        x = x+1
        current_pos = model[center_index].position
        model = translation(model, axis=1, surface=surf)
        if x > 9:
            break
    # Align in x-direction
    x = 0
    while not a/2 > (current_pos[0] - center_atom_position[0]) > -a/2:
        x = x+1
        current_pos = model[center_index].position
        model = translation(model, axis=0, surface=surf)
        if x > 9:
            break

    # Reapply constraints ince operation is complete
    if prev_const is not None:
        model.set_constraint(prev_const)

    return model
