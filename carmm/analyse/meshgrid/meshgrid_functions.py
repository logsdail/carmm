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
        MeshObject: MeshgridObject
            Object storing meshgrid and PBC conditions
    Returns:
        mesh_distances: Numpy meshgrid of distances from point [a_xx, a_yy, a_zz].
    """

    import numpy as np

    # Convert supplied point to fractional coordinates.

    x_dist = np.abs(a_xx - MeshObject.xx)
    y_dist = np.abs(a_yy - MeshObject.yy)
    z_dist = np.abs(a_zz - MeshObject.zz)

    if MeshObject.pbc[0]:
        x_dist = np.where(x_dist > 0.5 * MeshObject.x_max, x_dist - MeshObject.x_max, x_vec_dist)
    if MeshObject.pbc[1]:
        y_dist = np.where(y_dist > 0.5 * MeshObject.y_max, y_dist - MeshObject.y_max, y_vec_dist)
    if MeshObject.pbc[2]:
        z_dist = np.where(z_dist > 0.5 * MeshObject.z_max, z_dist - MeshObject.z_max, z_vec_dist)

    mesh_distances = np.sqrt(x_dist**2 + y_dist**2 + z_dist**2)
    #    frac_a = np.dot(np.array([a_xx, a_yy, a_zz]),MeshObject.inverse_cell_array)

#    x_vec_dist = np.abs(frac_a[0] - MeshObject.frac_xx)
#    y_vec_dist = np.abs(frac_a[1] - MeshObject.frac_yy)
#    z_vec_dist = np.abs(frac_a[2] - MeshObject.frac_zz)

#    x_dist = np.abs(a_xx - MeshObject.xx)
#    y_dist = np.abs(a_yy - MeshObject.yy)
#    z_dist = np.abs(a_zz - MeshObject.zz)

#    if MeshObject.pbc[0]:
#        x_dist = np.where(x_vec_dist > 0.5, x_vec_dist - 1., x_vec_dist)
#    if MeshObject.pbc[1]:
#        y_dist = np.where(y_vec_dist > 0.5, y_vec_dist - 1., y_vec_dist)
#    if MeshObject.pbc[2]:
#        z_dist = np.where(z_vec_dist > 0.5, z_vec_dist - 1., z_vec_dist)



    # Convert mesh distances back to cartesian coordinates.
    #stack_mesh = np.stack((x_dist, y_dist, z_dist), axis=-1)
    #stack_mesh = np.einsum('ji,abcj->iabc', MeshObject.Cell.array, stack_mesh)

    #mesh_distances = np.linalg.norm(stack_mesh, axis=0)

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
        MeshObject: MeshgridObject
            Object storing meshgrid and PBC conditions
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
        MeshObject: MeshgridObject
            Object storing meshgrid and PBC conditions
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

    mic_shift = find_mic((vec2 - vec1), cell=MeshObject.Cell.array)

    midpoint = np.linalg.norm(vec1+(mic_shift/2))

    return midpoint

def atom_mesh_build_mask(MeshObject, Atom):
    """
    Defines meshgrid of based on a radius set by atomic_radius. Uses 99.99 as a junk value.
    FUTURE DEVELOPMENT:
    - Replace junk value with None??
    Args:
        MeshObject: MeshgridObject
            Object storing meshgrid and PBC conditions
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
        mol_yy = np.where(new_distances < atom_radius, MeshObject.yy, mol_yy)
        mol_zz = np.where(new_distances < atom_radius, MeshObject.zz, mol_zz)

    return mol_xx, mol_yy, mol_zz
