def void_mesh_build_coarse_grid(atoms, void_min, void_max, ucell_obj, mol_xx, coarseness=1):
    # Defines the mesh points occupied by atoms.
    X0 = np.full(ucell.nx, fill_value=99.99)
    Y0 = np.full(ucell.ny, fill_value=99.99)
    Z0 = np.full(ucell.nz, fill_value=99.99)
    void_xx, void_yy, void_zz = np.meshgrid(X0, Y0, Z0, indexing='xy')
    void_centres = []

    # FIND UNOCCUPIED SITES
    probe_xx = np.where(mol_xx == 99.99, ucell.xx, 99.99)

    mol_xx_mask = (mol_xx != 99.99)

    continue_count = 0
    not_continue_count = 0

    for i in tqdm(range(0, ucell.nx, coarseness)):
        for j in range(0, ucell.ny, coarseness):
            for k in range(0, ucell.nz, coarseness):

                if probe_xx[i, j, k] == 99.99:
                    continue_count += 1
                    continue

                not_continue_count += 1
                a_xx, a_yy, a_zz = ucell_obj.xx[i, j, k], ucell_obj.yy[i, j, k], ucell_obj.zz[i, j, k]

                point_accepted = False
                for pbox_r in np.arange(void_min, void_max, 0.5):
                    probe_box_distance_matrix = distance_MIC_mesh(a_xx, a_yy, a_zz, ucell)

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