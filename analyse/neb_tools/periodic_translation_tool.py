def translation(model):
    '''
    Takes X arguments: filename/Atoms object,
    TODO: manipulate cell, rotate and change position of atoms
    consistent with symmetry - WITHOUT changing indices

    Parameters:
    model: Atoms object or string
        If string, e.g.: 'name.traj', a file of this name will be read
        to retrieve model.
    '''

    from ase.io import read
    import numpy as np
    import math
    from software.analyse.neb_tools.check_geometry import switch_indices

    def sort_by_xyz(model):
        '''
        Sorting indices by xyz coordinates
        Returns sorted index list
        '''
        from ase.build import sort
        sort(model)
        xyz = model.get_positions()
        sorted_xyz = sorted(xyz, key=lambda k: [k[2], k[1]])
        index_by_xyz = []
        current_index = -1

        # Find sequence in which atoms are now arranged by xyz
        for positions in xyz: # [x,y,z]
            current_index += 1
            count_iter = - 1  # Indexing in Atoms object starts from 0 not 1
            for sorted_positions in sorted_xyz:
                count_iter += 1
                if (sorted_positions == positions).all():
                    index_by_xyz = index_by_xyz + [[current_index, count_iter]]

        # Find the sequence to switch indices correctly
        index_by_xyz = sorted(index_by_xyz, key=lambda k:[k[1]])
        change_xyz = []

        for switch in index_by_xyz:
            print(switch[0])
            change_xyz += [switch[0]]

        return change_xyz

    if isinstance(model, str) is True:
        model = read(model)

    '''Section on variables'''
    # XY SHIFT
    # take 2.0 A worth of atoms to the other side. then all cordinates - 2.0A
    a = 3.914
    shift_dist = a * (2**(1/2)/2)
    # Take into account the  deg rotation, include tolerance for surf atoms
    shift_coordinate = shift_dist * math.cos(math.radians(60))
    #print(shift_coordinate)
    indices_to_move = []
    # Retrieve constraints from the model for later
    constraint = model._get_constraints()
    indices_old = sort_by_xyz(model)


    '''Section on moving atoms'''
    # align atoms perpendicular to x-axis
    model.rotate(30, 'z', rotate_cell=True)

    positions = model.get_positions()
    count_iter = 0
    for coordinate in positions:
        count_iter = count_iter + 1
        if coordinate[0] < shift_coordinate:
            # #print(coordinate[0])
            indices_to_move = indices_to_move + [count_iter - 1]

    model.rotate(-30, 'z', rotate_cell=True)

    # rpt cell, take from one side and add temp, exchange change positions
    temp_model = model.repeat((2, 1, 1))
    temp_model.set_constraint()
    model.set_constraint()
    # temp_model.set_calculator(model.get_calculator())

    temp_indices = np.array(indices_to_move) + len(model.get_tags())
    #print(indices_to_move)
    # Merge list of indices to switch
    index_pairs = zip(indices_to_move, temp_indices)

    for pairs in index_pairs:
        model[pairs[0]].position = temp_model[pairs[1]].position

    model.set_positions(
                    np.array(model.get_positions()) - np.array(
                                                        (shift_dist, 0, 0)))

    '''Section on index correction'''
    model.rotate(30, 'z', rotate_cell=True)
    indices_new = sort_by_xyz(model)
    index_pairs = zip(indices_new, indices_old)

    model = model[indices_new]

    # model.set_constraint(constraint)
    return model
