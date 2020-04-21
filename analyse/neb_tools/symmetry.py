'''This file is work in progress'''


def translation(model, a, axis=0, surface="111"):
    '''
    Performs a translaton of the model by manipulation of the unit cell,
    maintaining the optimised geometry and forces. After translation original
    Atoms object is permanently changed.

    TODO:
    - Is forces array shifting with the atoms in the Atoms object?
            Adjustment added in sort_by_xyz - test required to make sure it is
            performed twice.
    - FOR NOW requires surface to have tags for layers of atoms in Z-direction
    - functionality beyond FCC? or higher index?
    - provide a reasonable example

    Parameters:
    model: Atoms object or string
        If string, e.g.: 'name.traj', a file of this name will be read
        to retrieve model.
    a: float
        lattice parameter used in the model
    axis: integer
        Choice - 0, 1 representing axis x, y
    surface: string
        FCC surface - so far supports "111", "110", "100"
    '''

    from ase.io import read
    import numpy as np
    import math

    if isinstance(model, str) is True:
        model = read(model)

    # Retrieve constraints, calculator from the model for later
    constraint = model._get_constraints()
    prev_calc = model.get_calculator()

    '''Section on variables'''
    indices_to_move = []
    # TODO: Make lattice parameter independent based on closest neighbours
    shift_dist = a * (2**(1/2)/2)
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
                                (a, 0, 0)))
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

    # Retain calculator information and constraint
    if model.get_calculator() is not None:
        prev_calc.atoms = model
    model.set_calculator(prev_calc)
    model.set_constraint(constraint)

    return model


def sort_by_xyz(model, surface):
    ''' WORK IN PROGRESS
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


def mirror(model, center_index, plane="y", surf="111"):
    '''
    Function that returns a mirror image of the adsorbate in the x or y axis
    and shifts the surface atoms accordingly.

    Parameters:
    model: Atoms object
        periodic FCC surface model with an adsorbate. Tags required for surface
        layers.
    plane: string
        "x" or "y"
    center_index: integer
        Index of an atom with respect to which the image will be
        mirrored.
    surf: string
        Surface type. Supported "111", "110", "100" of the FCC lattice

    # TODO: turn into a function that can take either x or y and a specific
    # index of an adsorbate atom as the center
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

    no_atoms = len(model.get_tags())
    translate = model.positions[center_index][axis]

    for i in [atom.index for atom in model]:
        model.positions[i][axis] = (-model.positions[i][axis] + (2*translate))
    # align to position zero
    zero = (model[[atom.index for atom in model if atom.tag == 1][0]].position[axis])

    for i in reversed([atom.index for atom in model]):
        model.positions[i][axis] = (model.positions[i][axis] - zero)

    def borrow_positions(model, axis, surf):
        ''' TODO: detach and reuse elsewhere '''
        from software.analyse.neb_tools.symmetry import sort_by_xyz

        indices_to_move = []
        constraint = model._get_constraints()
        count_iter = 0
        # align axis
        if surf == "111":
            model.rotate(30, 'z', rotate_cell=True)
        positions = model.get_positions()

        # Mark atoms that need moving
        for coordinate in positions:
            count_iter = count_iter + 1
            # Include some tolerance for surface atoms
            if coordinate[axis] < -1.0:
               indices_to_move = indices_to_move + [count_iter - 1]
        # Return cell to original shape
        if surf == "111":
            model.rotate(-30, 'z', rotate_cell=True)

        if axis == 1:
            temp_model = model.repeat((1, 2, 1))
        elif axis == 0:
            temp_model = model.repeat((2, 1, 1))

        temp_model.set_constraint()
        model.set_constraint()
        for i in indices_to_move:
            model[i].position[axis] = temp_model[i+len(
                model.get_tags())].position[axis]

        # Add retrieved constraints, calculator
        model.set_constraint(constraint)
        model = sort_by_xyz(model, surf)
        return model

    from ase.visualize import view
    # TODO: Execute x times to ensure no axis coordinate below -0.50
    model = borrow_positions(model, axis, surf)
    model = borrow_positions(model, axis, surf)

    if model.get_calculator() is not None:
        prev_calc.atoms = model
        model.set_calculator(prev_calc)

    from software.analyse.neb_tools.symmetry import sort_by_xyz
    model = sort_by_xyz(model, surf)

    return model
