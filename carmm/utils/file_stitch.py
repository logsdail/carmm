def file_stitch(path, out_fname, data, lines_per_image=None, mode='1toA'):

    """
    A small function that stitches an array of atomwise data (partial charge, displacement, force, etc) into an xyz file
    containing geometry data. Allows for construction of extxyz from regular xyz plus extra required data.
    Args:
        path: Path to xyz geometry file (str)
        out_fname: Path to output file (str)
        data: Array of data to be appended (array)
        lines_per_image: Number of lines per image if appending to xyz trajectory. This is number of atoms + 2 (int)
        mode: determines the mode of the function. 1toA (1 to All) appends a dataset that is n_atoms long to every set
              of atoms in an image, giving all atoms the dame labels. AtoA (All to All) assumes that the dataset is
              continuous over n_atoms * n_images, labelling all atoms with all data. (str)

    """

    import shutil
    import fileinput

    shutil.copyfile(path, out_fname)

    if lines_per_image is None:
        # File is n_atoms(1) + title(1) + n_coords(n) lines long
        n_lines_image = data.size + 2
    else:
        n_lines_image = lines_per_image

    image = 0

    with fileinput.input(out_fname, inplace=True) as p:
        for line in p:
            if (p.lineno() - 1) % n_lines_image == 0:
                # p.lineno() - 1 is used as lines are numbered from 1 and arrays are indexed from 0
                # Update image number
                image += 1
            if (p.lineno() - 1) % n_lines_image == 0 or (p.lineno() - 1) % n_lines_image == 1:
                # skip line 0 and 1 for each image to preserve .xyz structure
                print(line, end='')

            else:
                # edit file in place to append data. Formatting is defined by location of newline
                if mode == '1toA':
                    # assign n_atoms data to n_atoms * n_images lines
                    data_str = str(data[:][((p.lineno()-1) % n_lines_image) - 2])
                elif mode == 'AtoA':
                    # assign n_atoms * n_images data to n_atoms * n_images lines
                    # 2 lines are added per image as the n_atoms and title lines for each individual image
                    data_str = str(data[:][(p.lineno() - 1) - (2*image)])
                if '[' in data_str:
                    data_str = data_str.replace("[","")
                    data_str = data_str.replace("]","")
                mod = line[0:-1] + '   ' + data_str + '\n'
                print(mod, end='')
