
# Contains a function to stitch extra data onto the end of an existing .xyz file to create an extxyz file

def file_stitch(path, out_fname, data, lines_per_image=None):

    """

    Args:
        path: Path to xyz geometry file (str)
        out_fname: Path to output file (str)
        data: Array of data to be appended (array)
        lines_per_image: Number of lines per image if appending to xyz trajectory (int)

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
                # Update image number
                image += 1
            if (p.lineno() - 1) % n_lines_image == 0 or (p.lineno() - 1) % n_lines_image == 1:
                # skip line 0 and 1 for each image to preserve .xyz structure
                print(line, end='')

            else:
                # edit file in place to append data. Formatting is defined by location of newline
                mod = line[0:-1] + '   ' + str(data[((p.lineno() - 1) % n_lines_image) - 2]) + '\n'
                print(mod, end='')
