def extract_mulliken_charge(fn, natoms):
    '''
    Function to extract and return the Mulliken charges from an FHI-aims output.
    A list of relative charges is returned, with +q meaning charg depletion and -q meaning charge accummulation

    Parameters:

    fn: string
        Filename from which the Mulliken data should be extracted
    natoms: int
        Number of atoms in the calculation

    '''

    with open(fn, 'r') as f:
        output = f.readlines()

    # Focus on just the Mulliken data
    mulliken_line = 0
    output_line = len(output) - 1
    while output_line and not mulliken_line:
        if "Starting Mulliken Analysis" in output[output_line]:
            mulliken_line = output_line
            # 8 lines from text label to start of data (spin-paired)
            mulliken_line += 8
        else:
            output_line -= 1

    mulliken_data = [q.split()[3] for q in output[mulliken_line:mulliken_line + natoms]]

    return mulliken_data