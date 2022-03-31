def distance_meshgrid2point(a_xx, a_yy, a_zz, unit_cell_object, mic):
    '''
    Function finds the distance between a point (defined in cartesian co-ordinates a_xx, a_yy, a_zz)
    relative to all points on a Numpy meshgrid (defined in the unit_cell_object).
    Args:
        a_xx, a_yy, a_zz: float
            Cartesian co-ordinates of point.
        unit_cell_object: Ucell object
            Object storing parameters of the meshgrid.
        MIC: logical
            Determines whether the minimum image convention should be used.
    Returns:
        mesh_distances: Numpy meshgrid of distances from point [a_xx, a_yy, a_zz].
    '''

    import numpy as np

    mesh_distances = np.sqrt(
        (a_xx - unit_cell_object.xx) ** 2 + (a_yy - unit_cell_object.yy) ** 2 + (a_zz - unit_cell_object.zz) ** 2)

    if (mic):
        for x in np.arange(-1, 2):
            for y in np.arange(-1, 2):
                for z in np.arange(-1, 2):
                    new_distances = np.sqrt((a_xx - unit_cell_object.xx + (unit_cell_object.dim[0] * x)) ** 2
                                            + (a_yy - (unit_cell_object.yy + (unit_cell_object.dim[1] * y))) ** 2
                                            + (a_zz - (unit_cell_object.zz + (unit_cell_object.dim[2] * z))) ** 2)

                    mesh_distances = np.where(new_distances < mesh_distances, new_distances, mesh_distances)

    return mesh_distances


def distance_point2point(x_1, y_1, z_1, x_2, y_2, z_2, mic, unit_cell):
    '''
    Function finds the distance between two points (defined in cartesian co-ordinates).
    Args:
        x_1, y_1, z_1: float
            Cartesian co-ordinates of point 1.
        x_2, y_2, z_2: float
            Cartesian co-ordinates of point 2.
        MIC: logical
            Determines whether the minimum image convention should be used.
    Returns:
        o_distance: Distance between points 1 and 2.
    '''

    import numpy as np

    o_distance = np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)

    if (mic):
        for x in np.arange(-1, 2):
            for y in np.arange(-1, 2):
                for z in np.arange(-1, 2):
                    new_distance = np.sqrt((x_1 - (x_2 + (unit_cell.dim[0] * x)) ** 2) +
                                           (y_1 - (y_2 + (unit_cell.dim[1] * y)) ** 2) +
                                           (z_1 - (z_2 + (unit_cell.dim[2] * z)) ** 2))
                    if o_distance > new_distance:
                        o_distance = new_distance

    return o_distance

def midpoint_MIC_points(x_1, y_1, z_1, x_2, y_2, z_2, MIC, dim):
    '''
    Function finds the distance between two points (defined in cartesian co-ordinates).
    Args:
        x_1, y_1, z_1: float
            Cartesian co-ordinates of point 1.
        x_2, y_2, z_2: float
            Cartesian co-ordinates of point 2.
        MIC: logical
            Determines whether the minimum image convention should be used.
    Returns:
        x_mid, y_mid, z_mid: Midpoint between point 1 and 2.
    '''

    import numpy as np

    o_distance = np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)
    x_mid, y_mid, z_mid = (x_1 + x_2) / 2, (y_1 + y_2) / 2, (z_1 + z_2) / 2

    if (MIC):
        for x in np.arange(-1, 2):
            for y in np.arange(-1, 2):
                for z in np.arange(-1, 2):
                    x_sft = dim[0] * x
                    y_sft = dim[1] * y
                    z_sft = dim[2] * z
                    new_distance = np.sqrt((x_1 - (x_2 + x_sft)) ** 2 + (y_1 - (y_2 + y_sft)) ** 2 + ((z_2 + z_sft)) ** 2)
                    if o_distance > new_distance:
                        o_distance = new_distance

                        x_mid, y_mid, z_mid = (x_1 + (x_2 + x_sft)) / 2, (y_1 + (y_2 + y_sft)) / 2, (
                                    z_1 + (z_2 + z_sft)) / 2
                        if not 0 <= x_mid <= dim[0]:
                            x_mid = x_mid - x_sft
                        if not 0 <= y_mid <= dim[1]:
                            y_mid = y_mid - y_sft
                        if not 0 <= z_mid <= dim[2]:
                            z_mid = z_mid - z_sft

    return x_mid, y_mid, z_mid

def atom_mesh_build_mask(ucell, atom, atom_radius, mic):
    '''
    Defines meshgrid of based on a radius set by atomic_radius. Uses 99.99 as a junk value.
    FUTURE DEVELOPMENT:
    - Atom by atom definitions of atom_radius.
    - Replace junk value with None??
    Args:
        ucell: unit_cell object
            Contains meshgrid and unit cell information
        atom: Atoms object
            Contains information (atomic symbols and positions)
        atom_radius: float
            Determines whether the minimum image convention should be used.
        mic: logical
            Minimum image convention for PBC on/off.
    Returns:
        mol_xx, mol_yy, mol_zz: Contains values of (x, y, z) co-ordinates for each point on axis occupied by
        atom van der Waals volume. Set to junk value otherwise.
    '''

    import numpy as np

    # Defines the mesh points occupied by atoms. Uses 99.99 as a trash value.
    X0 = np.full(ucell.nx, fill_value=99.99)
    Y0 = np.full(ucell.ny, fill_value=99.99)
    Z0 = np.full(ucell.nz, fill_value=99.99)

    mol_xx, mol_yy, mol_zz = np.meshgrid(X0, Y0, Z0, indexing='xy')

    for atom_co in range(len(atom.positions)):
        a_xx, a_yy, a_zz = atom.positions[atom_co][0], atom.positions[atom_co][1], atom.positions[atom_co][2]

        new_distances = distance_meshgrid2point(a_xx, a_yy, a_zz, ucell, mic)

        mol_xx = np.where(new_distances < atom_radius, ucell.xx, mol_xx)
        mol_yy = np.where(new_distances < atom_radius, ucell.yy, mol_yy)
        mol_zz = np.where(new_distances < atom_radius, ucell.zz, mol_zz)

    return mol_xx, mol_yy, mol_zz

def atom_mesh_build_mask_pbox(ucell, atom, atom_radius, mic):
    '''
    Defines meshgrid of based on a radius set by atomic_radius. Uses 99.99 as a junk value.
    FUTURE DEVELOPMENT:
    - Atom by atom definitions of atom_radius.
    - Replace junk value with None??
    Args:
        ucell: unit_cell object
            Contains meshgrid and unit cell information
        atom: Atoms object
            Contains information (atomic symbols and positions)
        atom_radius: float
            Determines whether the minimum image convention should be used.
        mic: logical
            Minimum image convention for PBC on/off.
    Returns:
        mol_xx, mol_yy, mol_zz: Contains values of (x, y, z) co-ordinates for each point on axis occupied by
        atom van der Waals volume. Set to junk value otherwise.
    '''

    import numpy as np

    # Defines the mesh points occupied by atoms. Uses 99.99 as a trash value.
    X0 = np.full(ucell.nx, fill_value=99.99)
    Y0 = np.full(ucell.ny, fill_value=99.99)
    Z0 = np.full(ucell.nz, fill_value=99.99)

    mol_xx, mol_yy, mol_zz = np.meshgrid(X0, Y0, Z0, indexing='xy')

    for atom_co in range(len(atom.positions)):

        a_xx, a_yy, a_zz = atom.positions[atom_co][0], atom.positions[atom_co][1], atom.positions[atom_co][2]
        act_boxes, act_boxes_mic = find_active_boxes(a_xx,a_yy,a_zz,ucell,atom_radius,mic)

        for b in range(len(act_boxes)):
            bx, by, bz = act_boxes[b][0],act_boxes[b][1],act_boxes[b][2]

            x_mic, y_mic, z_mic = act_boxes_mic[b][0], act_boxes_mic[b][1], act_boxes_mic[b][2]

            bxmin, bymin, bzmin = ucell.ind_list_x[bx][0], ucell.ind_list_y[by][0], ucell.ind_list_z[bz][0]
            bxmax, bymax, bzmax = ucell.ind_list_x[bx][1], ucell.ind_list_y[by][1], ucell.ind_list_z[bz][1]

            new_distances = np.sqrt((a_xx - (ucell.xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax] + (ucell.dim[0] * x_mic))) ** 2
                                    + (a_yy - (ucell.yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax] + (ucell.dim[1] * y_mic))) ** 2
                                    + (a_zz - (ucell.zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax] + (ucell.dim[2] * z_mic))) ** 2)

            mol_xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(new_distances < atom_radius,
                                                                     ucell.xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax],
                                                                     mol_xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax])
            mol_yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(new_distances < atom_radius,
                                                                     ucell.yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax],
                                                                     mol_yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax])
            mol_zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(new_distances < atom_radius,
                                                                     ucell.zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax],
                                                                     mol_zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax])

    return mol_xx, mol_yy, mol_zz

def find_active_boxes(x1, y1, z1, unit_cell, radius, mic):

    import numpy as np

    active_boxes = []
    active_boxes_mic = []

    half_diag=np.sqrt(unit_cell.part_dx**2+unit_cell.part_dy**2+unit_cell.part_dz**2)/2

    for ind in range(len(unit_cell.box_means)):
        print(unit_cell.box_means[ind][0])
        print(unit_cell.box_means[ind][1])
        print(unit_cell.box_means[ind][2])
        point2part_dist, x_mic, y_mic, z_mic = distance_point2point(x1, y1, z1, unit_cell.box_means[ind][0],
                                                                   unit_cell.box_means[ind][1],
                                                                   unit_cell.box_means[ind][2], mic, unit_cell)
        if (point2part_dist-radius) < half_diag:
            active_boxes.append([unit_cell.box_mean_indices[ind][0], unit_cell.box_mean_indices[ind][1],
                                 unit_cell.box_mean_indices[ind][2]])
            active_boxes_mic.append([x_mic, y_mic, z_mic])

    return active_boxes, active_boxes_mic