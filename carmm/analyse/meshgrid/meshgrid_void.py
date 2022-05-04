def void_find(meshobject, atoms, coarseness=1):
    """
    Defines the maximum spherical radius of probes in unoccupied space. Scans through the distance
    between the probe point and all atoms in the system, and defines the maximum radius of the probe
    as the smallest value of (probe_position - atom_position - vdw_radii).
    Coarseness factor reduces number of probe points by skipping over points in the grid.
    Args:
        meshobject: Mesh object
            Contains meshgrid and associated unit cell information.
        atoms: Atoms object
            Supplies coordinates and atomic information.
        coarseness: integer
            Skips over an integer number of points
    Return:
        void_centres: numpy array
            Atomic coordiantes of the void centres.
        void_radii: numpy array
            Radius of the void centres
    """

    import numpy as np
    from ase.data.vdw_alvarez import vdw_radii
    from ase.geometry import get_distances

    # Check atoms PBC and mesh PBC match.
    if meshobject.strict_mode:
        mol_mesh_pbc_check(mesh.pbc, atoms.pbc)

    void_centres = []
    void_radii = []

    atom_radii = vdw_radii[atoms.numbers]

    atom_radii = np.array([atom_radii]).flatten()

    for i in range(0, meshobject.nx, coarseness):
        for j in range(0, meshobject.ny, coarseness):
            for k in range(0, meshobject.nz, coarseness):

                a_xx, a_yy, a_zz = meshobject.xx[i, j, k], meshobject.yy[i, j, k], meshobject.zz[i, j, k]

                distances = get_distances(np.array([a_xx, a_yy, a_zz]), p2=atoms.positions,
                                          cell=meshobject.Cell, pbc=meshobject.pbc)

                distances = distances[1] - atom_radii

                min_distance = np.amin(distances)

                void_centres.append([a_xx, a_yy, a_zz])
                void_radii.append([min_distance])

    void_centres = np.array([void_centres])
    void_centres = void_centres[0]
    void_radii = np.array([void_radii]).flatten()

    return void_centres, void_radii


def void_build_mask(meshobject, void_centres, void_radii, min_void=1.0):
    """
    Defines meshgrid of voids given by list of void centres and their radii. Uses np.nan as a junk value.
    Args:
        meshobject: Mesh object.
            Meshgrid and unit cell information.
        void_centres: numpy array
            Coordinates of the void centres.
        void_radii: numpy array
            Radius of the void sphere.
        min_void: float
            Defines the minimum size of void to plot.
    Returns:
        void_xx, void_yy, void_zz: numpy array
            Contains values of (x, y, z) coordinates for each point on axis occupied by
            atom van der Waals volume. Set to junk value otherwise.
    """

    import numpy as np
    from carmm.analyse.meshgrid.meshgrid_functions import distance_meshgrid2point

    # Defines the mesh points occupied by atoms. Uses np.nan as a trash value.
    x0 = np.full(meshobject.nx, fill_value=np.nan)
    y0 = np.full(meshobject.ny, fill_value=np.nan)
    z0 = np.full(meshobject.nz, fill_value=np.nan)

    void_xx, void_yy, void_zz = np.meshgrid(x0, y0, z0, indexing='xy')

    for i in range(len(void_centres)):

        if void_radii[i] < min_void:
            continue

        a_xx, a_yy, a_zz = void_centres[i][0], void_centres[i][1], void_centres[i][2]

        new_distances = distance_meshgrid2point(a_xx, a_yy, a_zz, meshobject)

        void_xx = np.where(new_distances < void_radii[i], meshobject.xx, void_xx)
        void_yy = np.where(new_distances < void_radii[i], meshobject.yy, void_yy)
        void_zz = np.where(new_distances < void_radii[i], meshobject.zz, void_zz)

    return void_xx, void_yy, void_zz


def void_analysis(meshobject, void_centres, void_radii, void_xx):
    """
    Prints basic volumentric information of the void.
    Args:
        meshobject: Mesh object.
            Meshgrid and unit cell information.
        void_centres: numpy array
            Coordinates of the void centres.
        void_radii: numpy array
            Radius of the void sphere.
        void_xx: meshgrid numpy array
            x coordinates of void mesh grid.
    Returns:
        void_volume: float
            The total volume of the volume accessible to the probe (Ang**3).
    """

    import numpy as np

    largest = np.argmax(void_radii)

    print(f"Maximum void radius of {void_radii[largest]} at {void_centres[largest]}")

    unique, counts = np.unique(void_xx, return_counts=True)
    counter = counts[-1]

    volume = meshobject.Cell.volume

    unocc_sites = np.size(void_xx) - counter
    vox_volume = volume / (meshobject.nx * meshobject.ny * meshobject.nz)

    void_volume = unocc_sites * vox_volume

    print(f"Total volume of void: {void_volume} Ang**3")
    print(f"Total volume of unit cell: {volume} Ang**3")

    return void_volume
