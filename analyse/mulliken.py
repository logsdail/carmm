# TODO: Duplicate this for spin.
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

#TODO: Document this like crazy as it is super unclear at present
def store_k_point_data(previous_k_point, current_k_point, current_atom, previous_atom_added,
        all_k_points_energies, all_k_points_occupancies, all_k_points_mulliken, energies, occupancy, mulliken):
    '''

    :param previous_k_point:
    :param current_k_point:
    :param current_atom:
    :param previous_atom_added:
    :param all_k_points_energies:
    :param all_k_points_occupancies:
    :param all_k_points_mulliken:
    :param energies:
    :param occupancy:
    :param mulliken:
    :return:
    '''

    # Store k-point information when looking at 1st atom
    if current_atom == 1 or previous_atom_added == 1 or (current_atom == 2 and previous_atom_added == 0):
        if previous_k_point > len(all_k_points_energies):
            all_k_points_energies.append(energies)
            all_k_points_occupancies.append(occupancy)
        else:
            all_k_points_energies[previous_k_point-1] = all_k_points_energies[previous_k_point-1] + energies
            all_k_points_occupancies[previous_k_point-1] = all_k_points_occupancies[previous_k_point-1] + occupancy

    # Collect all k-point information for Mulliken. Array structure: all_k_points_mulliken[k_point][atom][mulliken]
    if previous_k_point > len(all_k_points_mulliken):
        all_k_points_mulliken.append([0 for x in range(0)])

    #print previous_k_point, current_k_point, previous_atom_added, current_atom
    if previous_k_point >= current_k_point and current_atom != previous_atom_added and previous_atom_added > 0:
        if previous_atom_added > len(all_k_points_mulliken[previous_k_point-1]):
            #print "Appending1: ", previous_k_point, previous_atom_added, len(mulliken), len(all_k_points_mulliken[previous_k_point-1])
            all_k_points_mulliken[previous_k_point-1].append(mulliken)
        else:
            #print "Inserting1: ", previous_k_point, previous_atom_added, len(mulliken), len(all_k_points_mulliken[previous_k_point-1])
            all_k_points_mulliken[previous_k_point-1][previous_atom_added-1] = [ old + new for old, new in zip(all_k_points_mulliken[previous_k_point-1][previous_atom_added-1], mulliken) ]
    else:
        if current_atom > len(all_k_points_mulliken[previous_k_point-1]):
            #print "Appending2: ", previous_k_point, previous_atom_added, len(mulliken), len(all_k_points_mulliken[previous_k_point-1])
            all_k_points_mulliken[previous_k_point-1].append(mulliken)
        else:
            #print "Inserting2: ", previous_k_point, previous_atom_added, len(mulliken), len(all_k_points_mulliken[previous_k_point-1])
            all_k_points_mulliken[previous_k_point-1][current_atom-1] = [ old + new for old, new in zip(all_k_points_mulliken[previous_k_point-1][current_atom-1], mulliken) ]

#TODO: document this like crazy
def parse_mulliken_file(input_lines):
    '''

    :param input_lines:
    :return:
    '''

    # Initialise variables
    mulliken = [0 for x in range(0)]
    all_k_points_energies = [0 for x in range(0)]
    all_k_points_occupancies = [0 for x in range(0)]
    all_k_points_weights = [0 for x in range(0)]
    all_k_points_mulliken = [0 for x in range(0)]
    current_atom = 0
    previous_atom_added = 0
    current_k_point = 0
    previous_k_point = 0

    for sentence in input_lines:
        split_sentence = sentence.split()
        if len(split_sentence) > 0:
            if split_sentence[0] == "Atom":
                current_atom = int(split_sentence[2].replace(":",""))
            elif split_sentence[0] == "k":
                previous_k_point = current_k_point #Just added - to be used to improve storage of all data for all k-points
                current_k_point = int(split_sentence[3].replace(":",""))
                # Store k-point weight. Compatible only with 180327 build onwards
                if len(split_sentence) > 9:
                    all_k_points_weights.append(float(split_sentence[10]))

                if previous_k_point > 0:
                    store_k_point_data(previous_k_point, current_k_point, current_atom, previous_atom_added,
                        all_k_points_energies, all_k_points_occupancies, all_k_points_mulliken, energies, occupancy, mulliken)
                    previous_atom_added = current_atom

                #Reset arrays
                energies = [0 for x in range(0)]
                occupancy = [0 for x in range(0)]
                mulliken = [[0 for x in range(0)] for x in range(6)]
            elif split_sentence[0] == "Spin" or split_sentence[0] == "#" or split_sentence[0] == "State":
                continue
            else:
                # Read in all data for this line and put in the correct places
                energies.append(float(split_sentence[1]))
                occupancy.append(float(split_sentence[2]))
                # This line is kind of redundant, as we can calculate it from angular contributions but it does save some recalculation
                mulliken[0].append(float(split_sentence[3]))
                mulliken[1].append(float(split_sentence[4])) #s
                if len(split_sentence) > 5:
                    mulliken[2].append(float(split_sentence[5])) #p
                    if len(split_sentence) > 6:
                        mulliken[3].append(float(split_sentence[6])) #d
                        if len(split_sentence) > 7:
                            mulliken[4].append(float(split_sentence[7])) #f
                            if len(split_sentence) > 8:
                                mulliken[5].append(float(split_sentence[8])) #g

    # Needed to ensure dataset is complete for last k-point
    store_k_point_data(current_k_point, current_k_point, current_atom, previous_atom_added,
           all_k_points_energies, all_k_points_occupancies, all_k_points_mulliken, energies, occupancy, mulliken)

    # In case k-point weights are not defined i.e. pre 180319 output
    if len(all_k_points_weights) == 0:
        all_k_points_weights = [1.0]*len(all_k_points_energies)

    return all_k_points_weights, all_k_points_energies, all_k_points_occupancies, all_k_points_mulliken