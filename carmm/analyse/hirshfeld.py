# Much recycled from mulliken.py, this should get Hirshfeld charges from aims.out

def extract_hirshfeld(fname, natoms, data, write=True, outname='hirshfeld.txt'):
    """

    Args:

        fname: Input file name, usually aims.out. str
        natoms: number of atoms to extract charge from. int
        data: type of output requested. str
              Can be 'charge', 'volume', 'volume f', 'dipole vector', 'dipole moment' or 'second'
        write: bool of whether to write the hirshfeld data to a new file called hirshfeld.txt

    Returns: hirsh. List of requested data

    """

    # Extracts data from an aims.out and writes it to a new file.
    import sys

    ids = {
        'charge':'Hirshfeld charge        :',
        'volume':'Hirshfeld volume        :',
        'volume f':'Free atom volume        :',
        'dipole vector':'Hirshfeld dipole vector :',
        'dipole moment':'Hirshfeld dipole moment :',
        'second':'Hirshfeld second moments:'
    }

    if data in ids:
        identifier = ids.get(data)

    else:
        print('Requested data not recognised')
        sys.exit()

        
    with open(fname, 'r') as f:
        output = f.readlines()

    hirshfeld_line = get_hirshfeld_line(output)

    if write is True:

        write_hirshfeld(output, natoms, hirshfeld_line, outname)

        hirsh = read_hirshfeld(outname, identifier)

    else:
        hirsh = read_hirshfeld(fname, identifier)

    return hirsh


def get_hirshfeld_line(text):

    # Gets the line from which to read Hirshfeld data

    hirshfeld_line = 0
    output_line = len(text) - 1
    while output_line and not hirshfeld_line:
        if "Performing Hirshfeld analysis of fragment charges and moments." in text[output_line]:
            hirshfeld_line = output_line

        else:
            output_line -= 1

    return hirshfeld_line


def write_hirshfeld(text, natoms, start, outname):

    # Internal function to write out hirshfeld data into a new file
    # Text: list of lines from file
    # natoms: no. of atoms
    # start: line within file to start from
    # outname: name of file to write

    with open(outname, 'w') as h:

        line = start + 2
        for n in range(natoms):
            count = 0
            while count <= 9:
                # print(output[line+count])
                h.write(text[line + count])
                count += 1

            # Hirshfeld output is 10 lines long
            line = line + 10


def read_hirshfeld(fname, identifier):
    """

    Args:
        fname: Input filename. str
        identifier: regex to pull data from file. str

    Returns: label. list of data

    """
    label = list()

    with open(fname, 'r') as hf:
        data_lines = hf.readlines()

        for line in range(len(data_lines)):
            if identifier in data_lines[line]:
                if identifier == 'Hirshfeld second moments:':
                    hirsh_label = data_lines[line][32:-1].strip() + ' '
                    hirsh_label += data_lines[line+1][32:-1].strip() + ' '
                    hirsh_label += data_lines[line+2][32:-1].strip() + ' '
                else:
                    hirsh_label = data_lines[line][32:-1].strip()

                if identifier == 'Hirshfeld dipole vector :' or identifier == 'Hirshfeld second moments:':
                    internal = hirsh_label.split()
                    for i in range(len(internal)):
                        internal[i] = float(internal[i])
                    label.append(internal)
                else:
                    label.append(float(hirsh_label))

    return label


def vmd_out(array, fname='vmd_chrgs.txt'):

    # Quick script to write charges to a file for a TCL script input to place those charges on atoms in VMD
    # array = numpy array
    # fname is str of filename to write out

    with open(fname, 'w') as v:
        for line in range(array.size):
            v.write(str(array[line]))
            v.write('\n')

