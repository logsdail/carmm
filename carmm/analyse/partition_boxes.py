

def distance_MIC_partbox_mesh(a_x, a_y, a_z, ucell_obj, bxmin, bxmax, bymin, bymax, bzmin, bzmax, x_mic, y_mic, z_mic):
    new_distances = np.sqrt(
        (a_x - (ucell_obj.xx[bxmin:bxmax, bymin:bymax, bzmin:bzmax] + (ucell_obj.dim[0] * x_mic))) ** 2
        + (a_y - (ucell_obj.yy[bxmin:bxmax, bymin:bymax, bzmin:bzmax] + (ucell_obj.dim[1] * y_mic))) ** 2
        + (a_z - (ucell_obj.zz[bxmin:bxmax, bymin:bymax, bzmin:bzmax] + (ucell_obj.dim[2] * z_mic))) ** 2)

    return new_distances

def partition_indices(partition_axis, npoints, ind_list):
    p_per_box = npoints / partition_axis

    for npart in range(partition_axis):
        min_partind_ax = int(npart * p_per_box)
        max_partind_ax = int(((npart + 1) * p_per_box))

        ind_list.append([min_partind_ax, max_partind_ax])

    return ind_list

def calculate_part_box_means(p_x, p_y, p_z, i_listx, i_listy, i_listz):
    box_means = []
    box_mean_indices = []
    for x in range(p_x):
        for y in range(p_y):
            for z in range(p_z):
                x_mean = np.mean(
                    xx[i_listx[x][0]:i_listx[x][1], i_listy[y][0]:i_listy[y][1], i_listz[z][0]:i_listz[z][1]])
                y_mean = np.mean(
                    yy[i_listx[x][0]:i_listx[x][1], i_listy[y][0]:i_listy[y][1], i_listz[z][0]:i_listz[z][1]])
                z_mean = np.mean(
                    zz[i_listx[x][0]:i_listx[x][1], i_listy[y][0]:i_listy[y][1], i_listz[z][0]:i_listz[z][1]])
                box_mean_indices.append([x, y, z])
                box_means.append([x_mean, y_mean, z_mean])

    return box_means, box_mean_indices


def find_active_boxes(x1, y1, z1, unit_cell, alpha):
    active_boxes = []
    active_boxes_mic = []
    for ind in range(len(unit_cell.box_means)):
        point2part_dist, x_mic, y_mic, z_mic = distance_MIC_points(x1, y1, z1, unit_cell.box_means[ind][0],
                                                                   unit_cell.box_means[ind][1],
                                                                   unit_cell.box_means[ind][2], unit_cell.dim)
        if point2part_dist < alpha:
            active_boxes.append([unit_cell.box_mean_indices[ind][0], unit_cell.box_mean_indices[ind][1],
                                 unit_cell.box_mean_indices[ind][2]])
            active_boxes_mic.append([x_mic, y_mic, z_mic])

    return active_boxes, active_boxes_mic
