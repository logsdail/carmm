'''This file is work in progress'''
# TODO: Enable lattice parameter input - important for alloys
import copy


def translation(model, axis=0, surface="111", m_m_dist=None):
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

    assert axis in [0,1], "axis must be 0 for x or 1 for y"
    # Avoid changes to the original Atoms object
    model = copy.deepcopy(model)


    # Retrieve constraints from the model for later
    constraint = model._get_constraints()

    '''Section on variables'''
    indices_to_move = []
    # Extraction of lattice parameter from bottom layers
    if m_m_dist is not None:
        shift_dist = m_m_dist
    else:
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
    tolerance = 1.

    if axis == 0:
        # TODO: get rid of the rotation, it is there because math is easier
        #       would require editing sort_y_tag function
        # align atoms perpendicular to x-axis
        if surface == "111":
            model.rotate(30, 'z', rotate_cell=True)

        positions = model.get_positions()
        count_iter = 0
        for coordinate in positions:
            count_iter += 1
            # Include tolerance
            if coordinate[0] < shift_x*tolerance:
                indices_to_move += [count_iter - 1]

        if surface == "111":
            model.rotate(-30, 'z', rotate_cell=True)

        # rpt cell, take from one side and add temp, exchange change positions
        temp_model = model.repeat((2, 1, 1))
        temp_model.set_constraint()
        model.set_constraint()
        temp_indices = np.array(indices_to_move) + len(model)
        # Merge list of indices to switch
        index_pairs = zip(indices_to_move, temp_indices)

        for pairs in index_pairs:
            model[pairs[0]].position = temp_model[pairs[1]].position

        if surface in ["111", "100"]:
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
            if coordinate[1] < shift_y*tolerance:
                indices_to_move = indices_to_move + [count_iter - 1]

        # rpt cell, take from one side and add temp, exchange change positions
        temp_model = model.repeat((1, 2, 1))
        temp_model.set_constraint()
        model.set_constraint()
        temp_indices = np.array(indices_to_move) + len(model)

        # Merge list of indices to switch
        index_pairs = zip(indices_to_move, temp_indices)

        for pairs in index_pairs:
            model[pairs[0]].position = temp_model[pairs[1]].position

        if surface == "111":
            model.set_positions(
                            np.array(model.get_positions()) - np.array(
                                (shift_x, shift_y, 0)))
        elif surface in ["100", "110"]:
            model.set_positions(
                            np.array(model.get_positions()) - np.array(
                                (0, shift_y, 0)))

    model = wrap_fcc(model, surface)

    # Retain previous constraints
    model.set_constraint(constraint)

    return model


def wrap_fcc(model, surface):
    '''

    Args:
        model:
        surface:

    Returns:

    '''

    model = model.copy()
    # Retrieve constraints from the model for later
    constraint = model._get_constraints()

    '''Section on index correction'''
    model = sort_by_xyz(model, surface)

    zero_index = 0
    #zero_index = get_zero_from_constrained_atoms(model)

    zero_x = copy.deepcopy((model[zero_index].position[0]))
    zero_y = copy.deepcopy((model[zero_index].position[1]))

    # Reset to zero
    for i in reversed([atom.index for atom in model]):
        model.positions[i][0] = (model.positions[i][0] - zero_x)
        model.positions[i][1] = (model.positions[i][1] - zero_y)

    model = sort_by_xyz(model, surface)

    # Remove rounding errors and set corner to (0,0,z)
    for i in reversed([atom.index for atom in model]):
        model.positions[i][0] = (model.positions[i][0] - model[0].i)
        model.positions[i][1] = (model.positions[i][1] - model[0].y)


    model.wrap()
    model = sort_by_xyz(model, surface)

    # Retain previous constraints
    model.set_constraint(constraint)

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
        Periodic surface model, so far supported are monometallic
        FCC "111", "110", "100"
        HCP "0001"
    surface: string
        Face centered cubic low index surfce - "111", "110" or "110"
    '''
    import numpy as np
    # Avoid changes to the original Atoms object
    import copy
    model = copy.deepcopy(model)

    # Sorting mechanism for Y tags
    def sort_y_tag(xyz, surface):
        import numpy as np
        # TODO: shortest_AB could be extracted from bottom fixed layer for
        #    consistency, but poses problems if layers are not the same
        if surface in ["111", "0001"]:
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

    return model


def mirror(model, center_index, plane="y", surf="111", m_m_dist=None):
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
        Surface type. Supported "110", "100" of the FCC lattice
    '''
    # Avoid changes to the original Atoms object
    import copy
    model = copy.deepcopy(model)

    assert plane in ["x", "y"], "Wrong axis chosen. Please choose 'x' or 'y'"
    planes = {"x":0, "y":1}
    axis = planes[plane]

    constraint = model._get_constraints()
    model.set_constraint()

    # Save initial center atom position
    # Hypothesis: center_atom_position is a shallow reference that does not
    #    survive other operations unchanged.
    # Solution: a deepcopy required to retain functionality? Confirmed.
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
    model = sort_by_xyz(model, surf)

    zero_index = 0

    zero_x = (model[zero_index].position[0])
    zero_y = (model[zero_index].position[1])

    for i in reversed([atom.index for atom in model]):
        model.positions[i][0] = (model.positions[i][0] - zero_x - model.get_cell_lengths_and_angles()[0])
        model.positions[i][1] = (model.positions[i][1] - zero_y - model.get_cell_lengths_and_angles()[1])

    model = sort_by_xyz(model, surf)
    model = wrap_fcc(model, surf)
    model.set_constraint(constraint)



    return model


def rotate_fcc(model, surf, center_index=0, m_m_dist=None):
    '''
    Rotate FCC cell with respect to an atom by allowed increments, ie.
    "111" - 120 degrees
    "100" - 90 degrees
    "110" - 180 degrees

    Parameters:
    model:in Atoms object
        FCC low index surface, (111), (110) or (100)
    center_index: int
        index of an atom as the center of the rotation operation
    surf: str
        "111", "110", "100" allowed
    '''
    import numpy as np

    # Avoid changes to the original Atoms object
    model = copy.deepcopy(model)

    if surf == "111":
        degrees = -120
    elif surf == "100":
        degrees = -90
    elif surf == "110":
        degrees = -180
    else:
        pass
        # TODO Error message

    # Remove constraints for rotation operation
    prev_const = model._get_constraints()
    model.set_constraint()
    # this variable does not survive other operations unchanged, hence deepcopy
    center_atom_position = copy.deepcopy(model[center_index].position)
    model.rotate(degrees, 'z', center=center_atom_position)

    model = wrap_fcc(model, surf)

    # Reapply constraints ince operation is complete
    if prev_const is not None:
        model.set_constraint(prev_const)


    return model


def sort_z(model, diff=1):
    '''
    Function assigning tags to atomic layers within periodic surface models
    for easier identification and compatibility with other functionality in CARMM.

    Parameters:
    model: Atoms object
        Periodic surface model
    diff: float or int
       Expected z-distance between atomic layers in Angstroms
    TODO: Incorporate this functionality into other parts of the code that rely
    on z-tags
    '''
    import numpy as np
    import copy

    # avoid manipulation of the original Atoms object
    model = copy.deepcopy(model)

    # extract atomic positions and sort them by the z-coordinate
    # sort the Atoms objects within Atoms based on z-coordinate
    z_sorted = sorted(model, key=lambda k: k.position[2])

    bin_counter = 0

    atomic_layers = [[]]
    # try to compare each consecutive number  - check if within +/- diff from each other
    for i in range(0, len(z_sorted)):
        if not i == 0:
            if -diff < z_sorted[i].position[2] - z_sorted[i-1].position[2] < diff:
                atomic_layers[bin_counter] += [z_sorted[i]]
            else:
                bin_counter += 1
                atomic_layers += [[z_sorted[i]]]
        else:
            atomic_layers[bin_counter] += [z_sorted[i]]


    # reverse order for setting tags, i.e. top layer (adsorbate) is zero, first layer is 1 etc.
    # go through slices
    # Ensure top layer is tag 1 as in ase.build functions
    tag_counter = len(atomic_layers)+1
    for i in range(0, len(atomic_layers)):
        tag_counter -= 1
        # for each atoms in slice assign tag to the corresponding atom in the original Atoms object
        for j in atomic_layers[i]:
            model[j.index].tag = tag_counter

    return model
