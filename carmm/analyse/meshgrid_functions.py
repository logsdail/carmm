def distance_meshgrid(a_xx, a_yy, a_zz, unit_cell_object, MIC):
    '''
    Function finds the distance between a point (defined in cartesian co-ordinates a_xx, a_yy, a_zz)
    relative to all points on a Numpy meshgrid (defined in the unit_cell_object). Defined to
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

    mesh_distances = np.sqrt(
        (a_xx - unit_cell_object.xx) ** 2 + (a_yy - unit_cell_object.yy) ** 2 + (a_zz - unit_cell_object.zz) ** 2)

    if (MIC):
        for x in np.arange(-1, 2):
            for y in np.arange(-1, 2):
                for z in np.arange(-1, 2):
                    new_distances = np.sqrt((a_xx - unit_cell_object.xx + (unit_cell_object.dim[0] * x)) ** 2
                                            + (a_yy - (unit_cell_object.yy + (unit_cell_object.dim[1] * y))) ** 2
                                            + (a_zz - (unit_cell_object.zz + (unit_cell_object.dim[2] * z))) ** 2)

                    mesh_distances = np.where(new_distances < old_distances, new_distances, old_distances)

    return mesh_distances


def distance_MIC_points(x_1, y_1, z_1, x_2, y_2, z_2, dim):
    '''
    Function finds the distance between a point (defined in cartesian co-ordinates a_xx, a_yy, a_zz)
    relative to all points on a Numpy meshgrid (defined in the unit_cell_object). Defined to
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

    o_distance = np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)
    x_out, y_out, z_out = 0, 0, 0

    for x in np.arange(-1, 2):
        for y in np.arange(-1, 2):
            for z in np.arange(-1, 2):
                new_distance = np.sqrt((x_1 - (x_2 + (dim[0] * x))) ** 2 + (y_1 - (y_2 + (dim[1] * y))) ** 2 + (
                            z_1 - (z_2 + (dim[2] * z))) ** 2)
                if o_distance > new_distance:
                    o_distance = new_distance
                    x_out, y_out, z_out = x, y, z

    return o_distance, x_out, y_out, z_out

def midpoint_MIC_points(x_1, y_1, z_1, x_2, y_2, z_2, dim):
    o_distance = np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)
    x_mid, y_mid, z_mid = (x_1 + x_2) / 2, (y_1 + y_2) / 2, (z_1 + z_2) / 2
    x_out, y_out, z_out = 0, 0, 0
    new_distance = 0

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

                    x_out, y_out, z_out = x, y, z

    return x_mid, y_mid, z_mid, x_out, y_out, z_out

def atom_mesh_build_mask(ucell, atom_obj, atom_radius):
    # Defines the mesh points occupied by atoms.
    X0 = np.full(ucell.nx, fill_value=99.99)
    Y0 = np.full(ucell.ny, fill_value=99.99)
    Z0 = np.full(ucell.nz, fill_value=99.99)

    mol_xx, mol_yy, mol_zz = np.meshgrid(X0, Y0, Z0, indexing='xy')

    for atom_co in range(len(atom_obj.coords)):
        a_xx, a_yy, a_zz = atom_obj.coords[atom_co][0], atom_obj.coords[atom_co][1], atom_obj.coords[atom_co][2]

        new_distances = distance_MIC_mesh(a_xx, a_yy, a_zz, ucell)

        mol_xx = np.where(new_distances < atom_radius, ucell.xx, mol_xx)
        mol_yy = np.where(new_distances < atom_radius, ucell.yy, mol_yy)
        mol_zz = np.where(new_distances < atom_radius, ucell.zz, mol_zz)

    return mol_xx, mol_yy, mol_zz
