def void_find(ucell_obj, atom, mic=True, coarseness=1):
    """
    Defines the void centres and radii of points across the meshgrid. Coarseness factor reduces number of probe points
    by skipping over points in the grid. Uses 99.99 as a junk value.
    Args:
        ucell_obj: unit_cell object
            Contains meshgrid and associated unit cell information.
        atom: Atoms object
            Supplies coordinates and atomic information.
        mic: logical
            Use minimum image convention.
        coarseness: integer
            Skips over an integer number of points
    Return:
        void_centres: numpy array
            Atomic co-ordiantes of the void centres.
        void_radii: numpy array
            Radius of the void centres
    """

    import numpy as np
    from ase.data.vdw_alvarez import vdw_radii
    from ase.data import atomic_numbers

    void_centres = []
    void_radii = []
    atom_radii = []

    for atom_co in range(len(atom.positions)):
        at_symbols = atom.get_chemical_symbols()
        at_n = atomic_numbers[at_symbols[atom_co]]
        atom_radius = vdw_radii[at_n]
        atom_radii.append(atom_radius)

    atom_radii = np.array([atom_radii]).flatten()

    for i in range(0, ucell_obj.nx, coarseness):
        for j in range(0, ucell_obj.ny, coarseness):
            for k in range(0, ucell_obj.nz, coarseness):

                a_xx, a_yy, a_zz = ucell_obj.xx[i, j, k], ucell_obj.yy[i, j, k], ucell_obj.zz[i, j, k]

                distances = np.abs(atom.positions - np.array([a_xx, a_yy, a_zz]))
                if mic:
                    distances = np.where(distances > 0.5 * ucell_obj.dim, distances - ucell_obj.dim, distances)

                distances = np.sqrt((distances ** 2).sum(axis=-1)) - atom_radii

                min_distance = np.amin(distances)

                void_centres.append([a_xx, a_yy, a_zz])
                void_radii.append([min_distance])

    void_centres = np.array([void_centres])
    void_centres = void_centres[0]
    void_radii = np.array([void_radii]).flatten()

    return void_centres, void_radii


def void_build_mask(ucell, void_centres, void_radii, mic, min_void=1.0):
    """
    Defines meshgrid of voids given by list of void centres and their radii. Uses 99.99 as a junk value.
    - Replace junk value with None??
    Args:
        ucell: unit_cell object
            Meshgrid and unit cell information
        void_centres: numpy array
            Co-ordinates of the void centres.
        void_radii: numpy array
            Radius of the void sphere.
        mic: logical
            Minimum image convention for PBC on/off.
        min_void: float
            Defines the minimum size of void to plot.
    Returns:
        void_xx, void_yy, void_zz: numpy array
            Contains values of (x, y, z) co-ordinates for each point on axis occupied by
            atom van der Waals volume. Set to junk value otherwise.
    """

    import numpy as np
    from carmm.analyse.meshgrid_functions import distance_meshgrid2point

    # Defines the mesh points occupied by atoms. Uses 99.99 as a trash value.
    X0 = np.full(ucell.nx, fill_value=99.99)
    Y0 = np.full(ucell.ny, fill_value=99.99)
    Z0 = np.full(ucell.nz, fill_value=99.99)

    void_xx, void_yy, void_zz = np.meshgrid(X0, Y0, Z0, indexing='xy')

    for i in range(len(void_centres)):

        if void_radii[i] < min_void:
            continue

        a_xx, a_yy, a_zz = void_centres[i][0], void_centres[i][1], void_centres[i][2]

        new_distances = distance_meshgrid2point(a_xx, a_yy, a_zz, ucell.xx, ucell.yy, ucell.zz, ucell.dim, mic)

        void_xx = np.where(new_distances < void_radii[i], ucell.xx, void_xx)
        void_yy = np.where(new_distances < void_radii[i], ucell.yy, void_yy)
        void_zz = np.where(new_distances < void_radii[i], ucell.zz, void_zz)

    return void_xx, void_yy, void_zz

def void_analysis(ucell, void_centres, void_radii, void_xx, void_yy, void_zz, mic):
    """
    Defines meshgrid of voids given by list of void centres and their radii. Uses 99.99 as a junk value.
    - Replace junk value with None??
    Args:
        ucell: unit_cell object
            Meshgrid and unit cell information
        void_centres: numpy array
            Co-ordinates of the void centres.
        void_radii: numpy array
            Radius of the void sphere.
        mic: logical
            Minimum image convention for PBC on/off.
        min_void: float
            Defines the minimum size of void to plot.
    """

    import numpy as np

    largest=np.argmax(void_radii)

    print(f"Maximum void radius of {void_radii[largest]} at {void_centres[largest]}")

    unique, counts = np.unique(void_xx, return_counts=True)
    counter=counts[-1]

    unocc_sites=np.size(void_xx)-counter
    vox_volume=ucell.dim[0]/ucell.nx*ucell.dim[1]/ucell.ny*ucell.dim[2]/ucell.nz
    total_volume=ucell.dim[0]*ucell.dim[1]*ucell.dim[2]

    void_volume=unocc_sites*vox_volume

    print(f"Total volume of void: {void_volume} Ang**3")
    print(f"Total volume of unit cell: {total_volume} Ang**3")
    print(f"Total vdw volume: {total_volume-void_volume} Ang**3")
