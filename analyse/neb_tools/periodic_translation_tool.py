'''This file is work in progress'''
def sort_by_xyz_old(model):
    ''' NEEDS TO BE SEQUENTIAL BY TAGS INSTEAD OF Z AXIS SORTING
    Sorting indices by xyz coordinates
    Returns sorted index list
    '''
    from ase.build import sort
    sort(model)
    xyz = model.get_positions()
    sorted_xyz = sorted(xyz, key=lambda k: [k[2], k[1], k[0]])
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

    return model[change_xyz]


def sort_by_xyz(model):
    ''' WORK IN PROGRESS
    Sorting indices by xyz coordinates
    Returns sorted index list
    '''
    import numpy as np
    # sort z direction by tags
    # reverse, need to start from bottom, ie. last tag
    tags = list(set(model.get_tags()))
    tags = sorted(tags, reverse=True)
    index_by_xyz = []
    xyz_all = model.get_positions()

    for tag in tags:
        index_by_tag = [atom.index for atom in model if atom.tag == tag]
        xyz = model[index_by_tag].get_positions()

        if not tag == 0:
            all_AB = model[index_by_tag].get_all_distances()
            shortest_ABs = []
            # Figure minimum bond lengths in this set of atoms
            for distance in all_AB:
                # Need to remove 0.0 from array before finding min value
                distance = np.amin(np.setdiff1d(distance,np.array(0.0)))
                shortest_ABs += [distance]
                # print(shortest_ABs)
            def sort_y_tag(xyz):
                y_tag = np.int(np.round((xyz[1])/(
                            np.sin(np.radians(60))*np.amin(shortest_ABs))))
                return y_tag

            sorted_xyz = sorted(
                xyz, key=lambda k: [sort_y_tag(k),
                                    k[0]])  # Y/ABmin - integer
        else:
            sorted_xyz = sorted(
                xyz, key=lambda k: [k[1], k[0]])

        # check all positions by tags and append a list
        for coordinate in sorted_xyz:
            test = (np.int((coordinate[1])/(
                np.sin(np.radians(60))*np.amin(shortest_ABs))))
            #print(coordinate[0])
            #print(index_by_xyz)
            for index in index_by_tag:
                if (coordinate == model[index].position).all():
                    index_by_xyz += [index]

    sort_test = np.argsort(index_by_xyz)

    print(index_by_xyz)
    model = model[index_by_xyz]

    return model

def translation(model):
    '''
    Takes X arguments: filename/Atoms object,
    TODO: manipulate cell, rotate and change position of atoms
    consistent with symmetry - WITHOUT changing indices
    For now - specific to fcc111

    Parameters:
    model: Atoms object or string
        If string, e.g.: 'name.traj', a file of this name will be read
        to retrieve model.
    TODO:
    a: float
        lattice parameter
    surface: string
        e.g. "fcc111", "fcc110"
    '''

    from ase.io import read
    import numpy as np
    import math

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
    # indices_old = sort_by_xyz(model)


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
    model = sort_by_xyz(model)

    return model
