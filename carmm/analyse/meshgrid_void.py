def void_mesh_build(void_min, void_max, ucell_obj, mol_xx, mic, coarseness=1):
    '''
    Defines the void (unoccupied regions) of a unit cell on a numpy meshgrid. Scans over probe points on the
    underlying meshgrid by increasing the radius of the grid until the van der Waals' volume of a molecule is
    encountered. Coarseness factor reduces number of probe points by skipping over points in the grid.
    Uses 99.99 as a junk value.
    Args:
        void_min: float
            Minimum radius scanned from for probe sphere.
        void_max: float
            Maximum radius scanned from for probe sphere.
        ucell_obj: unit_cell object
            Contains meshgrid and associated unit cell information.
        mol_xx: numpy array [nx,ny,nz]
            Grid to determine whether molecule is present at a point.
        mic: logical
            Use miniumum image convention.
        coarseness: integer
            Skips over an integer number of points
    Return:
        void_xx, void_yy, void_zz: Meshgrid of (x,y,z) occupied by the probe spheres.
    '''

    import numpy as np
    from carmm.analyse.meshgrid_functions import distance_meshgrid2point

    # Defines the mesh points occupied by atoms.
    X0 = np.full(ucell_obj.nx, fill_value=99.99)
    Y0 = np.full(ucell_obj.ny, fill_value=99.99)
    Z0 = np.full(ucell_obj.nz, fill_value=99.99)
    void_xx, void_yy, void_zz = np.meshgrid(X0, Y0, Z0, indexing='xy')
    void_centres = []

    # FIND UNOCCUPIED SITES
    probe_xx = np.where(mol_xx == 99.99, ucell_obj.xx, 99.99)

    mol_xx_mask = (mol_xx != 99.99)

    for i in range(0, ucell_obj.nx, coarseness):
        for j in range(0, ucell_obj.ny, coarseness):
            for k in range(0, ucell_obj.nz, coarseness):

                if probe_xx[i, j, k] == 99.99:
                    continue

                a_xx, a_yy, a_zz = ucell_obj.xx[i, j, k], ucell_obj.yy[i, j, k], ucell_obj.zz[i, j, k]

                point_accepted = False
                for pbox_r in np.arange(void_min, void_max, 0.5):
                    probe_box_distance_matrix = distance_meshgrid2point(a_xx, a_yy, a_zz, ucell_obj, mic)

                    ping = probe_box_distance_matrix < pbox_r

                    point_denied = np.any(ping & mol_xx_mask)

                    if (not point_denied):
                        point_accepted = True

                    if point_denied and (not point_accepted):
                        #print(f"Probe denied at: {a_xx,a_yy,a_zz} with radius: {pbox_r}")
                        break

                    if point_denied and point_accepted:
                        #print(f"Probe accepted at: {a_xx,a_yy,a_zz} with radius: {pbox_r-0.5}")
                        void_centres.append([a_xx, a_yy, a_zz])

                        ping = probe_box_distance_matrix < (pbox_r - 0.5)

                        void_xx = np.where(ping, ucell_obj.xx, void_xx)
                        void_yy = np.where(ping, ucell_obj.yy, void_yy)
                        void_zz = np.where(ping, ucell_obj.zz, void_zz)
                        break

    return void_xx, void_yy, void_zz
