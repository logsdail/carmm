# Much recycled from mulliken.py, this should get Hirshfeld charges from aims.out

def extract_hirshfeld(fname, natoms, data):
    """

    Args:
        fname: Input file name, usually aims.out. str
        natoms: number of atoms to extract charge from. int
        data: type of output requested. str
              Can be 'charge', 'volume', 'volume f', 'dipole vector', 'dipole moment' or 'second'

    Returns: label. List of requested data

    """

    # Extracts data from an aims.out and writes it to a new file.

    if data == 'charge':
        identifier = 'Hirshfeld charge        :'
    elif data == 'volume':
        identifier = 'Hirshfeld volume        :'
    elif data == 'volume f':
        identifier = 'Free atom volume        :'
    elif data == 'dipole vector':
        identifier = 'Hirshfeld dipole vector :'
    elif data == 'dipole moment':
        identifier = 'Hirshfeld dipole moment :'
    elif data == 'second':
        identifier = 'Hirshfeld second moments:'
    else:
        print('Requested data not recognised')

    with open(fname, 'r') as f:
        output = f.readlines()

        hirshfeld_line = get_hirshfeld_line(output)

    with open('hirshfeld.txt', 'w') as h:

        line = hirshfeld_line + 2
        for n in range(natoms):
            count = 0
            while count <= 9:
                h.write(output[line + count])
                count += 1

            # Hirshfeld output is 10 lines long
            line = line + 10

    label = list()

    with open(fname, 'r') as hf:
        data_lines = hf.readlines()

        for line in range(len(data_lines)):
            if identifier in data_lines[line]:
                hirsh_label = data_lines[line][32:-1].strip()
                label.append(hirsh_label)

    return label


def get_hirshfeld_line(text):

    # Gets the line from which to read Hirshfeld data

    hirshfeld_line = 0
    output_line = len(output) - 1
    while output_line and not hirshfeld_line:
        if "Performing Hirshfeld analysis of fragment charges and moments." in output[output_line]:
            hirshfeld_line = output_line

        else:
            output_line -= 1

    return hirshfeld_line

