def void_mesh_build(void_min, void_max, resol, ucell_obj, mol_xx, mic, coarseness=1):
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
                for pbox_r in np.arange(void_min, void_max, resol):
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

                        ping = probe_box_distance_matrix < (pbox_r - resol)

                        void_xx = np.where(ping, ucell_obj.xx, void_xx)
                        void_yy = np.where(ping, ucell_obj.yy, void_yy)
                        void_zz = np.where(ping, ucell_obj.zz, void_zz)
                        break

    return void_xx, void_yy, void_zz

def void_find_simple(void_min, void_max, resol, ucell_obj, atomco, atom_rad, mic, coarseness=1):
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
    from carmm.analyse.meshgrid_functions import distance_point2point

    void_centres = []
    void_radii = []

    for i in range(0, ucell_obj.nx, coarseness):
        for j in range(0, ucell_obj.ny, coarseness):
            for k in range(0, ucell_obj.nz, coarseness):

                a_xx, a_yy, a_zz = ucell_obj.xx[i, j, k], ucell_obj.yy[i, j, k], ucell_obj.zz[i, j, k]

                distances=[]
                for atom in atomco:

                    dist,xmic,ymic,zmic=distance_point2point(a_xx,a_yy,a_zz,
                                                             atom.position[0],atom.position[1],atom.position[2],
                                                             mic,ucell_obj.dim)
                    distances.append([dist])
                    point_accepted = False
                    point_denied = False

                distances=np.asarray(distances)
                for pbox_r in np.arange(void_min, void_max, resol):
                    if pbox_r<(distances+atom_rad):
                        point_accepted = True
                    else:
                        point_denied = False

                    if not(point_denied) and (point_accepted):
                        continue

                    if point_denied and (not point_accepted):
                        #print(f"Probe denied at: {a_xx,a_yy,a_zz} with radius: {pbox_r}")
                        break

                    if point_denied and point_accepted:
                        #print(f"Probe accepted at: {a_xx,a_yy,a_zz} with radius: {pbox_r-0.5}")
                        void_centres.append([a_xx, a_yy, a_zz])
                        void_radii.append([pbox_r-resol])
                        break

    void_centres=np.asarray(void_centres)
    void_radii=np.asarray(void_radii)

    return void_centres, void_radii

def void_mesh_build_pbox(void_min, void_max, resol, ucell, mol_xx, mic, coarseness=1):
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

    from carmm.analyse.meshgrid_functions import find_active_boxes
    import numpy as np

    # Defines the mesh points occupied by atoms.
    X0 = np.full(ucell.nx, fill_value=99.99)
    Y0 = np.full(ucell.ny, fill_value=99.99)
    Z0 = np.full(ucell.nz, fill_value=99.99)
    void_xx, void_yy, void_zz = np.meshgrid(X0, Y0, Z0, indexing='xy')
    test_xx, test_yy, test_zz = np.meshgrid(X0, Y0, Z0, indexing='xy')
    void_centres = []

    # FIND UNOCCUPIED SITES
    probe_xx = np.where(mol_xx == 99.99, ucell.xx, 99.99)

    mol_xx_mask = (mol_xx != 99.99)

    for i in range(0, ucell.nx, coarseness):
        for j in range(0, ucell.ny, coarseness):
            for k in range(0, ucell.nz, coarseness):

                if probe_xx[i, j, k] == 99.99:
                    continue

                a_xx, a_yy, a_zz = ucell.xx[i, j, k], ucell.yy[i, j, k], ucell.zz[i, j, k]

                point_accepted = False

                act_boxes, act_boxes_mic = find_active_boxes(a_xx, a_yy, a_zz, ucell, pbox_r, mic)

                bx, by, bz = act_boxes[b][0], act_boxes[b][1], act_boxes[b][2]

                x_mic, y_mic, z_mic = act_boxes_mic[b][0], act_boxes_mic[b][1], act_boxes_mic[b][2]

                bxmin, bymin, bzmin = ucell.ind_list_x[bx][0], ucell.ind_list_y[by][0], ucell.ind_list_z[bz][0]
                bxmax, bymax, bzmax = ucell.ind_list_x[bx][1], ucell.ind_list_y[by][1], ucell.ind_list_z[bz][1]

                for b in range(len(act_boxes)):

                    for pbox_r in np.arange(void_min, void_max, resol):

                        probe_box_distance_matrix = np.sqrt(
                            (a_xx + (ucell.dim[0] * x_mic) - ucell.xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax]) ** 2
                            + (a_yy + (ucell.dim[0] * x_mic) - ucell.yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax]) ** 2
                            + (a_zz + (ucell.dim[0] * x_mic) - ucell.zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax]) ** 2)

                        ping = probe_box_distance_matrix < pbox_r

                        point_denied = np.any(ping & mol_xx_mask[bxmin:bxmax, bymin:bymax, bzmin:bzmax])

                        test_xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(ping,
                                                                                  ucell.xx[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax],
                                                                                  void_xx[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax])
                        test_yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(ping,
                                                                                  ucell.yy[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax],
                                                                                  void_yy[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax])
                        test_zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(ping,
                                                                                  ucell.zz[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax],
                                                                                  void_zz[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax])

                    if (not point_denied):
                        void_xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(ping,
                                                                                  ucell.xx[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax],
                                                                                  void_xx[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax])
                        void_yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(ping,
                                                                                  ucell.yy[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax],
                                                                                  void_yy[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax])
                        void_zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax] = np.where(ping,
                                                                                  ucell.zz[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax],
                                                                                  void_zz[bxmin:bxmax, bymin:bymax,
                                                                                  bzmin:bzmax])
                        point_accepted = True

                    if point_denied and (not point_accepted):
                        #print(f"Probe denied at: {a_xx,a_yy,a_zz} with radius: {pbox_r}")
                        break

                    if point_denied and point_accepted:
                        #print(f"Probe accepted at: {a_xx,a_yy,a_zz} with radius: {pbox_r-0.5}")
                        void_centres.append([a_xx, a_yy, a_zz])

                        act_boxes, act_boxes_mic = find_active_boxes(a_xx, a_yy, a_zz, ucell, pbox_r-resol, mic)

                        for b in range(len(act_boxes)):
                            bx, by, bz = act_boxes[b][0], act_boxes[b][1], act_boxes[b][2]

                            x_mic, y_mic, z_mic = act_boxes_mic[b][0], act_boxes_mic[b][1], act_boxes_mic[b][2]

                            bxmin, bymin, bzmin = ucell.ind_list_x[bx][0], ucell.ind_list_y[by][0], \
                                                  ucell.ind_list_z[bz][0]
                            bxmax, bymax, bzmax = ucell.ind_list_x[bx][1], ucell.ind_list_y[by][1], \
                                                  ucell.ind_list_z[bz][1]

                            probe_box_distance_matrix = np.sqrt(
                                (a_xx + (ucell.dim[0] * x_mic) - ucell.xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax]) ** 2
                                + (a_yy + (ucell.dim[0] * x_mic) - ucell.yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax]) ** 2
                                + (a_zz + (ucell.dim[0] * x_mic) - ucell.zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax]) ** 2)

                            ping = probe_box_distance_matrix < (pbox_r - resol)

                        break

    return void_xx, void_yy, void_zz
