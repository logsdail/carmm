def void_find(MeshObject, Atom, coarseness=1):
    """
    Defines the void centres and radii of points across the meshgrid. Coarseness factor reduces number of probe points
    by skipping over points in the grid. Uses 99.99 as a junk value.
    Args:
        MeshObject: MeshgridObject object
            Contains meshgrid and associated unit cell information.
        Atom: Atoms object
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

    void_centres = []
    void_radii = []

    atom_radii = vdw_radii[Atom.numbers]

    atom_radii = np.array([atom_radii]).flatten()

    for i in range(0, MeshObject.nx, coarseness):
        for j in range(0, MeshObject.ny, coarseness):
            for k in range(0, MeshObject.nz, coarseness):

                a_xx, a_yy, a_zz = MeshObject.xx[i, j, k], MeshObject.yy[i, j, k], MeshObject.zz[i, j, k]

                distances = get_distances(np.array([a_xx, a_yy, a_zz]), p2=Atom.positions,
                                       cell=MeshObject.Cell, pbc=MeshObject.pbc)

                distances = distances - atom_radii

                min_distance = np.amin(distances)

                void_centres.append([a_xx, a_yy, a_zz])
                void_radii.append([min_distance])

    void_centres = np.array([void_centres])
    void_centres = void_centres[0]
    void_radii = np.array([void_radii]).flatten()

    return void_centres, void_radii


def void_build_mask(MeshObject, void_centres, void_radii, min_void=1.0):
    """
    Defines meshgrid of voids given by list of void centres and their radii. Uses 99.99 as a junk value.
    - Future change - replace junk value with None??
    Args:
        Ucell: unit_cell object
            Meshgrid and unit cell information
        void_centres: numpy array
            Coordinates of the void centres.
        void_radii: numpy array
            Radius of the void sphere.
        mic: logical
            Minimum image convention for PBC on/off.
        min_void: float
            Defines the minimum size of void to plot.
    Returns:
        void_xx, void_yy, void_zz: numpy array
            Contains values of (x, y, z) coordinates for each point on axis occupied by
            atom van der Waals volume. Set to junk value otherwise.
    """

    import numpy as np
    from carmm.analyse.meshgrid.meshgrid_functions import distance_meshgrid2point

    # Defines the mesh points occupied by atoms. Uses 99.99 as a trash value.
    x0 = np.full(MeshObject.nx, fill_value=99.99)
    y0 = np.full(MeshObject.ny, fill_value=99.99)
    z0 = np.full(MeshObject.nz, fill_value=99.99)

    void_xx, void_yy, void_zz = np.meshgrid(x0, y0, z0, indexing='xy')

    for i in range(len(void_centres)):

        if void_radii[i] < min_void:
            continue

        a_xx, a_yy, a_zz = void_centres[i][0], void_centres[i][1], void_centres[i][2]

        new_distances = distance_meshgrid2point(a_xx, a_yy, a_zz, MeshObject)

        void_xx = np.where(new_distances < void_radii[i], MeshObject.xx, void_xx)
        void_yy = np.where(new_distances < void_radii[i], MeshObject.yy, void_yy)
        void_zz = np.where(new_distances < void_radii[i], MeshObject.zz, void_zz)

    return void_xx, void_yy, void_zz


def void_analysis(MeshObject, void_centres, void_radii, void_xx):
    """
    Defines meshgrid of voids given by list of void centres and their radii. Uses 99.99 as a junk value.
    - Replace junk value with None??
    Args:
        Ucell: unit_cell object
            Meshgrid and unit cell information
        void_centres: numpy array
            Coordinates of the void centres.
        void_radii: numpy array
            Radius of the void sphere.
        void_xx: meshgrid numpy array
            x coordinates of void mesh grid.
    """

    import numpy as np

    largest = np.argmax(void_radii)

    print(f"Maximum void radius of {void_radii[largest]} at {void_centres[largest]}")

    unique, counts = np.unique(void_xx, return_counts=True)
    counter = counts[-1]

    xdim, ydim, zdim = MeshObject.array[0][0], MeshObject.array[1][1], MeshObject.array[2][2]

    unocc_sites = np.size(void_xx) - counter
    vox_volume = xdim / MeshObject.nx * ydim / MeshObject.ny * zdim / MeshObject.nz
    total_volume = xdim * ydim * zdim

    void_volume = unocc_sites * vox_volume

    print(f"Total volume of void: {void_volume} Ang**3")
    print(f"Total volume of unit cell: {total_volume} Ang**3")
    print(f"Total vdw volume: {total_volume - void_volume} Ang**3")
