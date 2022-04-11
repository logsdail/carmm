def distance_meshgrid2point(a_xx, a_yy, a_zz, xx, yy, zz, dim, mic):
    """
    Function finds the distance between a point (defined in cartesian co-ordinates a_xx, a_yy, a_zz)
    relative to all points on a Numpy meshgrid (defined in the unit_cell_object).
    Args:
        a_xx: float
            x coordinate of input point.
        a_yy: float
            y coordinate of input point.
        a_zz: float
            z coordinate of input point.
        xx: numpy array
            Meshgrid of x values on {nx,ny,nz} array.
        yy: numpy array
            Meshgrid of y values on {nx,ny,nz} array.
        zz: numpy array
            Meshgrid of z values on {nx,ny,nz} array.
        dim: numpy array
            Dimensions along each unit cell axis.
        mic: logical
            Minimum image convention On/Off.
    Returns:
        mesh_distances: Numpy meshgrid of distances from point [a_xx, a_yy, a_zz].
    """

    import numpy as np

    x_dist = np.abs(a_xx - xx)
    y_dist = np.abs(a_yy - yy)
    z_dist = np.abs(a_zz - zz)

    if mic:
        x_dist = np.where(x_dist > 0.5 * dim[0], x_dist - dim[0], x_dist)
        y_dist = np.where(y_dist > 0.5 * dim[1], y_dist - dim[1], y_dist)
        z_dist = np.where(z_dist > 0.5 * dim[2], z_dist - dim[2], z_dist)

    mesh_distances = np.sqrt(x_dist ** 2 + y_dist ** 2 + z_dist ** 2)

    return mesh_distances


def distance_point2point(x_1, y_1, z_1, x_2, y_2, z_2, mic, dim):
    """
    Function finds the distance between two points (defined in cartesian co-ordinates).
    Args:
        x_1: float
            x coordinate of point 1.
        y_1: float
            y coordinate of point 1.
        z_1: float
            z coordinate of point 1.
        x_2: float
            x coordinate of point 2.
        y_2: float
            y coordinate of point 2.
        z_2: float2
            z coordinate of point 2.
        mic: logical
            Determines whether the minimum image convention should be used.
        dim: numpy array
            Dimensions of unit cell in x, y, z.
    Returns:
        o_distance: Distance between points 1 and 2.
    """

    import numpy as np

    o_distance = np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)

    if mic:
        for x in np.arange(-1, 2):
            for y in np.arange(-1, 2):
                for z in np.arange(-1, 2):
                    new_distance = np.sqrt((x_1 - (x_2 + (dim[0] * x)) ** 2) +
                                           (y_1 - (y_2 + (dim[1] * y)) ** 2) +
                                           (z_1 - (z_2 + (dim[2] * z)) ** 2))
                    if o_distance > new_distance:
                        o_distance = new_distance

    return o_distance


def midpoint_points(x_1, y_1, z_1, x_2, y_2, z_2, mic, dim):
    """
    Function finds the distance between two points (defined in cartesian co-ordinates).
    Args:
        x_1: float
            x coordinate of point 1.
        y_1: float
            y coordinate of point 1.
        z_1: float
            z coordinate of point 1.
        x_2: float
            x coordinate of point 2.
        y_2: float
            y coordinate of point 2.
        z_2: float2
            z coordinate of point 2.
        mic: logical
            Determines whether the minimum image convention should be used.
        dim: numpy array
            x, y, z dimensions of the unit cell.
    Returns:
        x_mid: float
            x coordinate of midpoint between points 1 and 2.
        y_mid: float
            y coordinate of midpoint between points 1 and 2.
        z_mid: float
            z coordinate of midpoint between points 1 and 2.
    """

    import numpy as np

    o_distance = np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)
    x_mid, y_mid, z_mid = (x_1 + x_2) / 2, (y_1 + y_2) / 2, (z_1 + z_2) / 2

    if mic:
        for x in np.arange(-1, 2):
            for y in np.arange(-1, 2):
                for z in np.arange(-1, 2):
                    x_sft = dim[0] * x
                    y_sft = dim[1] * y
                    z_sft = dim[2] * z
                    new_distance = np.sqrt(
                        (x_1 - (x_2 + x_sft)) ** 2 + (y_1 - (y_2 + y_sft)) ** 2 + (z_2 + z_sft) ** 2)
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


def atom_mesh_build_mask(Ucell, atom, mic):
    """
    Defines meshgrid of based on a radius set by atomic_radius. Uses 99.99 as a junk value.
    FUTURE DEVELOPMENT:
    - Replace junk value with None??
    Args:
        Ucell: unit_cell object
            Contains meshgrid and unit cell information
        atom: Atoms object
            Contains information (atomic symbols and positions)
        mic: logical
            Minimum image convention for PBC on/off.
    Returns:
        mol_xx: 3D numpy array
            Value of x coordinate on {nx,ny,nz} grid.
        mol_yy: 3D numpy array
            Value of y coordinate on {nx,ny,nz} grid.
        mol_yy: 3D numpy array
            Value of z coordinate on {nx,ny,nz} grid.
    """

    import numpy as np
    from ase.data import atomic_numbers
    from ase.data.vdw_alvarez import vdw_radii

    # Defines the mesh points occupied by atoms. Uses 99.99 as a trash value.
    x0 = np.full(Ucell.nx, fill_value=99.99)
    y0 = np.full(Ucell.ny, fill_value=99.99)
    z0 = np.full(Ucell.nz, fill_value=99.99)

    mol_xx, mol_yy, mol_zz = np.meshgrid(x0, y0, z0, indexing='xy')

    for atom_co in range(len(atom.positions)):

        at_symbols = atom.get_chemical_symbols()
        at_n = atomic_numbers[at_symbols[atom_co]]
        atom_radius = vdw_radii[at_n]

        a_xx, a_yy, a_zz = atom.positions[atom_co][0], atom.positions[atom_co][1], atom.positions[atom_co][2]

        new_distances = distance_meshgrid2point(a_xx, a_yy, a_zz, Ucell.xx, Ucell.yy, Ucell.zz, Ucell.dim, mic)

        mol_xx = np.where(new_distances < atom_radius, Ucell.xx, mol_xx)
        mol_yy = np.where(new_distances < atom_radius, Ucell.yy, mol_yy)
        mol_zz = np.where(new_distances < atom_radius, Ucell.zz, mol_zz)

    return mol_xx, mol_yy, mol_zz
