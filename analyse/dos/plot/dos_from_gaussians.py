#!/usr/bin/env python3

def store_k_point_data(previous_k_point, current_k_point, current_atom, previous_atom_added, 
        all_k_points_energies, all_k_points_occupancies, all_k_points_mulliken, energies, occupancy, mulliken):

    # Store k-point information when looking at 1st atom
    if current_atom == 1 or previous_atom_added == 1:
        if previous_k_point > len(all_k_points_energies):
            all_k_points_energies.append(energies)
            all_k_points_occupancies.append(occupancy)
        else:
            all_k_points_energies[previous_k_point-1] = all_k_points_energies[previous_k_point-1] + energies
            all_k_points_occupancies[previous_k_point-1] = all_k_points_occupancies[previous_k_point-1] + occupancy

    #Collect all k-point information for Mulliken. Array structure: all_k_points_mulliken[k_point][atom][mulliken]
    if previous_k_point > len(all_k_points_mulliken):
        all_k_points_mulliken.append([0 for x in range(0)])

#    print previous_k_point, current_k_point, previous_atom_added, current_atom

    if previous_k_point >= current_k_point and current_atom != previous_atom_added and previous_atom_added > 0:
        all_k_points_mulliken[previous_k_point-1][previous_atom_added-1] = [ old + new for old, new in zip(all_k_points_mulliken[previous_k_point-1][previous_atom_added-1], mulliken) ]
    else:
        if current_atom > len(all_k_points_mulliken[previous_k_point-1]):
            all_k_points_mulliken[previous_k_point-1].append(mulliken)
        else:
            all_k_points_mulliken[previous_k_point-1][current_atom-1] = [ old + new for old, new in zip(all_k_points_mulliken[previous_k_point-1][current_atom-1], mulliken) ]

#    print len(all_k_points_energies[previous_k_point-1]), len(all_k_points_occupancies[previous_k_point-1])

def parse_data_file(input_file):

    # Initialise variables
    energies  = [[0 for x in range(0)]*1]
    occupancy = [[0 for x in range(0)]*1]
    
    for sentence in input_lines:
        #print sentence
        split_sentence = sentence.split()
        #print split_sentence

        energies[0].append(float(split_sentence[0]))
        occupancy[0].append(float(split_sentence[1]))

    #dummy variables
    k_weights = [1.0]*1
    atomic_mulliken = None

    #return energies[0], occupancy[0]
    return k_weights, energies, occupancy, atomic_mulliken

def parse_mulliken_file(input_file):

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
                #Storing all k-points will work if I just get the referencing correct
                #note that for instance the ordering for options in this if statement are currently incorrect
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
                mulliken = [[0 for x in range(0)] for x in range(5)] 
            elif split_sentence[0] == "Spin" or split_sentence[0] == "#" or split_sentence[0] == "State":
                continue 
            else:
                energies.append(float(split_sentence[1]))
                occupancy.append(float(split_sentence[2]))
                # This is kind of redundant, as we can calculate it from angular contributions
                mulliken[0].append(float(split_sentence[3])) 
                mulliken[1].append(float(split_sentence[4])) #s
                if len(split_sentence) > 5:
                    mulliken[2].append(float(split_sentence[5])) #p
                    if len(split_sentence) > 6:
                        mulliken[3].append(float(split_sentence[6])) #d
                        if len(split_sentence) > 7:
                            mulliken[4].append(float(split_sentence[7])) #f

    # Needed to ensure dataset is complete
    store_k_point_data(current_k_point, current_k_point, current_atom, previous_atom_added, 
           all_k_points_energies, all_k_points_occupancies, all_k_points_mulliken, energies, occupancy, mulliken)

#    for l in range(len(all_k_points_mulliken)):
#        for m in range(len(all_k_points_mulliken[l])):
#            print l, m, len(all_k_points_mulliken[l][m]), len(all_k_points_mulliken[l][m][0])

    # In case k-point weights are not defined
    if len(all_k_points_weights) == 0:
        all_k_points_weights = [1.0]*len(all_k_points_energies)

    return all_k_points_weights, all_k_points_energies, all_k_points_occupancies, all_k_points_mulliken

def calculate_homo_and_lumo_for_all_k_points(energies,occupancies):

    # Initialise variables
    homo = -999999.9
    lumo = 999999.9

    for kpt in range(len(occupancies)):
        temp_homo, temp_lumo = calculate_homo_and_lumo(energies[kpt],occupancies[kpt])
        print("k-point: ", kpt+1,"; HOMO: ", temp_homo, "; LUMO: ", temp_lumo, "; HOMO-LUMO: ", temp_homo-temp_lumo)
        if temp_homo > homo: 
            homo = temp_homo
        if lumo > temp_lumo:
            lumo = temp_lumo
  
    print("Overall: HOMO: ", homo, "; LUMO: ", lumo, "; HOMO-LUMO: ", homo-lumo)

    return homo, lumo

def calculate_homo_and_lumo(energies,occupancy):

    # Initialise variables
    homo = -999999.9
    lumo = 999999.9

    for j in range(len(occupancy)-1):
        #print occupancy[j]
        # Check if the occupancy of this state is more than 1, and less than one in the next state
        if occupancy[j] > 0.5 and occupancy[j+1] < 0.5:
            # Check if we already have a HOMO value that is higher, e.g. spin-polarised systems
            if energies[j] > homo:
                homo = energies[j]
            if energies[j+1] < lumo:
                lumo = energies[j+1]
#            print "XC: ", xc,"; HOMO: ", homo, "; LUMO: ", lumo, "HOMO-LUMO: ", homo-lumo

    return homo, lumo

def calculate_plot_values(min_point, max_point, separate_spin, energies, weights):

    import numpy as np
    from scipy.stats import norm

    points = 1000
    plot_values = [[0.0]*points]*2 #spin-up and spin-down
    x = np.linspace(min_point,max_point,points)

    #last_energy = energies[0][0]
    #spin_counter = 1

    #print xc, spin_counter, last_energy
    for kpt in range(len(energies)):
        #print weights[kpt]
        last_energy = energies[kpt][0]
        spin_counter = 1
        for e in range(len(energies[kpt])):
            energy = energies[kpt][e]
            if separate_spin and last_energy > energy:
                spin_counter += 1
                #print xc, spin_counter, last_energy, energy

            last_energy = energy
            if energy > min_point and energy < max_point:
                variance = 0.02
                sigma = np.sqrt(variance)
                plot_values[spin_counter-1] += weights[kpt]*norm.pdf(x,energy,sigma)

    #print plot_values

    return x, plot_values, spin_counter

def normalise_plot_values(plot_values, spin_counter, weight=1.0):

    for j in range(spin_counter):
        # Normalise
        max_value = max(plot_values[j])/weight
        #print j, max_value
        #plot_values /= max_value
        #No longer works so manually normalising over the array
        prefactor = 1
        if j > 0:
            prefactor = -1
        plot_values[j][:] = [prefactor*y/max_value for y in plot_values[j]]

###############################

import matplotlib.pyplot as plt
import sys
import pylab
import argparse

pylab.rcParams.update({'font.size': 24})
pylab.rcParams.update({'mathtext.fontset': 'stix'})
#plt.xlabel(r'$\epsilon$ (eV)')
#plt.ylabel('Density of States (1/eV)')

parser = argparse.ArgumentParser()
parser.add_argument('input_files', action='store',nargs="+",type=str,help="input filenames containing data tuples of energy and occupation")
parser.add_argument('--interpolate', action='store_true',help="interpolate under graph")
parser.add_argument('--lines', action='store_true',help="draw line graph")
parser.add_argument('--align_to_homo', action='store_true',help="align x-axis values relative to homo")
parser.add_argument('--split_updown', action='store_true',help="split spin-up and spin-down contributions")
parser.add_argument('--letters', action='store_true',help="give letters labelling each graph")
parser.add_argument('--xc_functionals', action='store',nargs="+",type=str,help="list of xc functional labels to iterate over for the filenames")
parser.add_argument('--shades', action='store_true',help="colour in shades instead of different colours")
parser.add_argument('--atoms', action='store',nargs="+",type=int,help="supply atom indices to include in atomic projected analysis")
parser.add_argument('--angular_momenta', action='store',type=str,help="supply angular momenta (from spdf) for projected analysis")
parser.add_argument('--gamma_only', action='store_true', help="only plot the DOS for the gamma-point of the system")
args =  parser.parse_args()

labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']
colors = ['red','blue','green','yellow','orange','indigo','violet']
if args.shades:
    #colors = ['LightCoral','red','darkred','brown']
    colors = ['pink','violet','MediumOrchid','purple']
interpolate_colors = ['Gold', 'MediumOrchid','LightCoral', 'LightGreen','LightPink','LightSkyBlue']
line_types = ['solid', 'dashed', 'dashdot', 'dotted']

ymax = 0
plt.figure(1)

if not args.xc_functionals:
    xc_functionals = [""]
else:
    xc_functionals = args.xc_functionals

k = -1
for xc in xc_functionals:
    input_files = [ i.replace(xc_functionals[0],xc) for i in args.input_files ]
    k = k + 1

    for i in range(len(input_files)):

        file = input_files[i]
        with open(file, 'r') as file_input:
            input_lines = file_input.readlines()

        if file[-12:] == "Mulliken.out":
            k_weights, energies, occupancies, all_atoms_total_mulliken = parse_mulliken_file(input_lines)
        else:
            k_weights, energies, occupancies, all_atoms_total_mulliken = parse_data_file(input_lines)

        # Dirty hack if the file doesn't contain any results
        if len(energies) == 0:
            break

        # Shrink arrays if only using the gamma-point
        if args.gamma_only:
            del k_weights[1:]
            del energies[1:]
            del occupancies[1:]
            if file[-12:] == "Mulliken.out":
                del all_atoms_total_mulliken[1:]

        print("Calculating HOMO/LUMO for: ", input_files[i])

        homo, lumo = calculate_homo_and_lumo_for_all_k_points(energies,occupancies)

        # Align everything to the homo?
        if args.align_to_homo:
            for kpt in range(len(energies)):
                energies[kpt][:] = [y-homo for y in energies[kpt]] 
            homo = 0.0
            plt.xlabel(r'$\epsilon - \epsilon_{HOMO}$ (eV)')

        #Subplots
        ax1 = plt.subplot(len(input_files),1,1)
        if i > 0:
            ax1 = plt.subplot(len(input_files),1,1+i,sharex=ax1)
        # Hack to set the y label where I want it
        ax1.set_yticklabels([])
        if i == len(input_files)/2:
            plt.ylabel('Density of States (1/eV)')

        min_point = homo - 20 
        max_point = homo + 20

        x, total_plot_values, spin_counter = calculate_plot_values(min_point,max_point,args.split_updown,energies,k_weights)

        normalise_plot_values(total_plot_values, spin_counter)

        print("Plotting graph ", i+1, " as ", labels[i], " in black : ", input_files[i])

        # Fill the background
        if args.interpolate and not (args.atoms or args.angular_momenta):
            for j in range(spin_counter):
                plt.fill_between(x, 0, total_plot_values[j],lw=0, facecolor=colors[i], label=input_files[i], interpolate=True, alpha=0.5)

        #Add on Mulliken analysis
        if all_atoms_total_mulliken and (args.atoms or args.angular_momenta):

            angular_momenta = {'a':0, 's':1, 'p':2, 'd':3, 'f':4}
            total_mulliken = [0 for x in range(0)]
            previous_total_mulliken = [0 for x in range(0)]  
            # Work out normalising constant
            mask = [0 for x in range(0)]
            total_electrons = [0 for x in range(0)]
            #total_electrons_alt = [0 for x in range(0)]
            for kpt in range(len(energies)):
                # Sometimes there are a different number of occupied states at each k-point. Is a hard cutoff of 0.5 electrons necessary?
                mask.append([(occ > 0.5) for occ in occupancies[kpt]])
                total_electrons.append(sum(mask[kpt]))
                #total_electrons_alt.append(sum(occupancies[kpt]))
                #print kpt, total_electrons[kpt], total_electrons_alt[kpt]
            #print len(energies)*max(total_electrons), sum(total_electrons), (sum(total_electrons_alt))

            if args.atoms:
                atoms_list = args.atoms
            else:
                atoms_list = range(1, len(all_atoms_total_mulliken[0])+1)

            if args.angular_momenta:
                angular_list = args.angular_momenta
            else:
                angular_list = "a"

            for angular in angular_list:
                #Store this in case we are interpolating
                if len(total_mulliken) > 0:
                    previous_total_mulliken = total_mulliken[:][:]
                 
                if angular not in angular_momenta:
                    print("Unknown option for angular momentum: ", angular)
                    print("Please choose from spdf only")
                    exit(1)

                # Effectively all we are doing here is normalising the values against the occupancy as a fraction of total electrons
                normalising_constant = 0.0
                for atom in atoms_list:
                    # Check if the specific angular momenta data we need actually exists for this atom by sampling Gamma-point
                    if len(all_atoms_total_mulliken[0][atom-1][angular_momenta[angular]]) > 0:
                        electrons_on_this_specific_atom = 0.0
                        for kpt in range(len(energies)):
                            electrons_on_this_specific_atom += sum([all_atoms_total_mulliken[kpt][atom-1][angular_momenta[angular]][occ] for occ in range(len(occupancies[kpt])) if mask[kpt][occ]])
                        # I don't really like the denominator here, but it seems to work
                        normalising_constant += electrons_on_this_specific_atom/(len(energies)*max(total_electrons))

                x, plot_values_mulliken, spin_counter = calculate_plot_values(min_point,max_point,args.split_updown,energies,k_weights)
                normalise_plot_values(plot_values_mulliken, spin_counter, normalising_constant)

                if len(total_mulliken) > 0:
                    total_mulliken = [ old + new for old, new in zip(total_mulliken, plot_values_mulliken) ]
                else:
                    total_mulliken = plot_values_mulliken[:][:]
                
                if len(previous_total_mulliken) == 0:
                    previous_total_mulliken = [[0] * len(tm) for tm in total_mulliken] 
 
                # Check the angular contribution is necessary to plot
                total_change = 0.0
                for j in range(spin_counter):
                    #print angular, sum(total_mulliken[j]), sum(previous_total_mulliken[j])
                    total_change = abs(sum(total_mulliken[j]) - sum(previous_total_mulliken[j]))
 
                    # Arbitrary threshold of half an electron:
                    if total_change > 0.5:
                        print("Plotting ", angular, " angular momenta in ", colors[angular_momenta[angular]-1], " for atoms ", atoms_list)
                        plt.plot(x,total_mulliken[j], lw=2, color=colors[angular_momenta[angular]-1], label=angular, ls=line_types[k])
                        if args.interpolate:
                            plt.fill_between(x, previous_total_mulliken[j], total_mulliken[j], lw=0, facecolor=colors[angular_momenta[angular]-1], label=input_files[i], interpolate=True)

        # Put this at the end so it covers everything else
        for j in range(spin_counter):
            plt.plot(x,total_plot_values[j],lw=2, color='black', label=input_files[i], ls=line_types[k])

        # Work to rescale axes
        # ycurrent = max(plot_values)+0.1*max(plot_values)
        # The above normalises to 1, so the below should always hold
        #ycurrent = 1.1
        #if ycurrent > ymax:
        ymax = 1.1 #ycurrent
        ymin = 0
        if args.split_updown:
            ymin = -ymax
        plt.ylim(ymin, ymax)
        plt.xlim(min_point+10, max_point-10)

        # HOMO
        plt.axvline(x=homo, ymin=-100, ymax=100, color='black', lw=2, ls=line_types[k]) # MFI

        # 0
        if args.split_updown:
            plt.axhline(y=0, xmin=-100, xmax=100, color='black', lw=2)

        # Label, hacked
        if args.letters:
            plt.text(max_point-16, 0.95, labels[i], verticalalignment='top')

#plt.legend(loc='upper right')

# Make sure this label is on the bottom graph
#ax1 = plt.subplot(len(input_files),1,len(input_files),sharex=ax1)
# Hack to set the y label where I want it
#ax1.set_yticklabels([])
plt.xlabel(r'$\epsilon$ (eV)')
#plt.ylabel('Density of States (1/eV)')
plt.show()

