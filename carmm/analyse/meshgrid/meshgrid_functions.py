def distance_meshgrid2point(a_xx, a_yy, a_zz, MeshObject):
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
    from ase.geometry import find_mic

    mesh_positions = np.stack((MeshObject.xx, MeshObject.yy, MeshObject.zz), axis=-1)
    mesh_positions = np.reshape(mesh_positions, ((MeshObject.nx * MeshObject.ny * MeshObject.nz), 3))

    distance_matrix = find_mic(np.array([a_xx,a_yy,a_zz]) - mesh_positions, MeshObject.Cell, MeshObject.pbc)

    mesh_distances = np.reshape(distance_matrix, (MeshObject.nx, MeshObject.ny, MeshObject.nz))

    return mesh_distances

def distance_point2point(x_1, y_1, z_1, x_2, y_2, z_2, MeshObject):
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

    from ase.geometry import get_distances

    o_distance = get_distances([x_1, y_1, z_1], p2=[x_2, y_2, z_2],
                                cell=MeshObject.Cell, pbc=MeshObject.pbc)

    return o_distance

def midpoint_points(x_1, y_1, z_1, x_2, y_2, z_2, MeshObject):
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
    from ase.geometry import find_mic

    vec1 = np.array([x_1, y_1, z_1])
    vec2 = np.array([x_2, y_2, z_2])

    mic_shift = find_mic((vec2 - vec1), cell=MeshObject.Cell.cell)

    midpoint = vec1+(mic_shift/2)

    return midpoint

def atom_mesh_build_mask(MeshObject, Atom):
    """
    Defines meshgrid of based on a radius set by atomic_radius. Uses 99.99 as a junk value.
    FUTURE DEVELOPMENT:
    - Replace junk value with None??
    Args:
        Ucell: unit_cell object
            Contains meshgrid and unit cell information
        Atom: Atoms object
            Contains information (atomic symbols and positions)
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
    x0 = np.full(MeshObject.nx, fill_value=99.99)
    y0 = np.full(MeshObject.ny, fill_value=99.99)
    z0 = np.full(MeshObject.nz, fill_value=99.99)

    mol_xx, mol_yy, mol_zz = np.meshgrid(x0, y0, z0, indexing='xy')

    for atom_co in range(len(Atom.positions)):

        at_symbols = Atom.get_chemical_symbols()
        at_n = atomic_numbers[at_symbols[atom_co]]
        atom_radius = vdw_radii[at_n]

        a_xx, a_yy, a_zz = Atom.positions[atom_co][0], Atom.positions[atom_co][1], Atom.positions[atom_co][2]

        new_distances = distance_meshgrid2point(a_xx, a_yy, a_zz, MeshObject)

        mol_xx = np.where(new_distances < atom_radius, MeshObject.xx, mol_xx)
        mol_yy = np.where(new_distances < atom_radius, MeshObject.xx, mol_yy)
        mol_zz = np.where(new_distances < atom_radius, MeshObject.xx, mol_zz)

    return mol_xx, mol_yy, mol_zz
