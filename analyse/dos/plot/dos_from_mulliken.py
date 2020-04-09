#!/usr/bin/env python3

# AJL, March 2018
# Script to take Mulliken.out and plot a density of states. Requires Matplotlib.
# For options, run as ./dos_from_mulliken.py -h
# For presentation effects (colours, font size, etc.) you will have to dig deeper below, roughly around lines 170 - 200
# Nov 2018: - updated print statements for Python3 compatbility

#### Analyse data files and work out our HOMO/LUMO or Fermi Level, depending on your preference

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
        # Check if the occupancy of this state is more than 1, and less than one in the next state
        if occupancy[j] > 0.5 and occupancy[j+1] < 0.5:
            # Check if we already have a HOMO value that is higher, e.g. spin-polarised systems
            if energies[j] > homo:
                homo = energies[j]
            if energies[j+1] < lumo:
                lumo = energies[j+1]

    return homo, lumo

#### Take energies and plot them as a DOS over the energy range of interest

def calculate_plot_values(min_point, max_point, separate_spin, energies, weights):

    import numpy as np
    from scipy.stats import norm

    points = 1000
    plot_values = [[0.0]*points]*2 #spin-up and spin-down
    x = np.linspace(min_point,max_point,points)

    for kpt in range(len(energies)):
        last_energy = energies[kpt][0]
        spin_counter = 1
        for e in range(len(energies[kpt])):
            energy = energies[kpt][e]
            if separate_spin and last_energy > energy:
                spin_counter += 1

            last_energy = energy
            if energy > min_point and energy < max_point:
                variance = 0.02
                sigma = np.sqrt(variance)
                plot_values[spin_counter-1] += weights[kpt]*norm.pdf(x,energy,sigma)

    return x, plot_values, spin_counter

def normalise_plot_values(plot_values, spin_counter, weight=1.0):

    for j in range(spin_counter):
        # Normalise
        max_value = max(plot_values[j])/weight
        #plot_values /= max_value
        #No longer works so manually normalising over the array
        prefactor = 1
        if j > 0:
            prefactor = -1
        plot_values[j][:] = [prefactor*y/max_value for y in plot_values[j]]

##############################

def exit_error(error_message):
    print("ERROR: "+error_message)
    exit(1)

###############################

#### Actual script

import matplotlib.pyplot as plt
import sys
import pylab
import argparse
# Mulliken output parser
from software.analyse.mulliken import parse_mulliken_file

pylab.rcParams.update({'font.size': 24})
pylab.rcParams.update({'mathtext.fontset': 'stix'})

parser = argparse.ArgumentParser()
parser.add_argument('input_files', action='store',nargs="+",type=str,help="input filenames containing data tuples of energy and occupation")
parser.add_argument('--align_to_homo', action='store_true',help="align x-axis values relative to homo")
parser.add_argument('--angular', action='store',type=str,help="supply angular momenta (in form e.g. 'spdfg') for projected analysis")
parser.add_argument('--atoms', action='store',nargs="+",type=int,help="supply atom indices to include in atomic projected analysis")
parser.add_argument('--gamma', action='store_true', help="only plot the DOS for the gamma-point of the system")
parser.add_argument('--interpolate', action='store_true',help="interpolate under graph")
parser.add_argument('--letters', action='store_true',help="give letters labelling each graph")
parser.add_argument('--polarised', action='store_true',help="split spin-up and spin-down contributions")
args =  parser.parse_args()

labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']
colors = ['red','blue','green','yellow','orange','indigo','violet']
line_types = ['solid', 'dashed', 'dashdot', 'dotted']
# This used to be a soft-coded variable that changed with filename, but for this script I've just disabled it
# If you want different line types, just change this to match the index from the list above that you want
k = 0 

ymax = 0
plt.figure(1)

for i in range(len(args.input_files)):

    file = args.input_files[i]
    with open(file, 'r') as file_input:
        input_lines = file_input.readlines()

    # Read Mulliken file
    k_weights, energies, occupancies, all_atoms_total_mulliken = parse_mulliken_file(input_lines)

    # Dirty hack if the file doesn't contain any results
    if len(energies) == 0:
        exit_error("No energies read. Error with input file?")

    # Shrink arrays if only using the gamma-point
    if args.gamma:
        del k_weights[1:]
        del energies[1:]
        del occupancies[1:]
        del all_atoms_total_mulliken[1:]

    print("Calculating HOMO/LUMO for: ", file) 

    homo, lumo = calculate_homo_and_lumo_for_all_k_points(energies,occupancies)

    # Align everything to the homo
    if args.align_to_homo:
        for kpt in range(len(energies)):
            energies[kpt][:] = [y-homo for y in energies[kpt]] 
        homo = 0.0

    # Subplots
    ax1 = plt.subplot(len(args.input_files),1,1)
    if i > 0:
        ax1 = plt.subplot(len(args.input_files),1,1+i,sharex=ax1)
    # Hack to set the y label where I want it
    ax1.set_yticklabels([])
    if i == len(args.input_files)/2:
        plt.ylabel('Density of States (1/eV)')

    # Arbitrary scale of interest for plot
    min_point = homo - 20 
    max_point = homo + 20

    # Plot values for DOS and normalise
    x, total_plot_values, spin_counter = calculate_plot_values(min_point,max_point,args.polarised,energies,k_weights)
    normalise_plot_values(total_plot_values, spin_counter)

    print("Plotting graph ", i+1, " as ", labels[i], " in black : ", file) 

    if args.interpolate and not (args.atoms or args.angular):
        # Fill the background
        for j in range(spin_counter):
            plt.fill_between(x, 0, total_plot_values[j], lw=0, facecolor=colors[i], label=file, interpolate=True, alpha=0.5)

    #Add on Mulliken analysis
    if all_atoms_total_mulliken and (args.atoms or args.angular):

        angular_momenta = {'a':0, 's':1, 'p':2, 'd':3, 'f':4, 'g':5}
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

        # Get atoms for projected DOS
        if args.atoms:
            atoms_list = args.atoms
            if max(atoms_list) > len(all_atoms_total_mulliken[0]):
                exit_error("Atom "+str(max(atoms_list))+" does not exist. Please check your atoms indices")
            elif min(atoms_list) < 1:
                exit_error("Atom "+str(min(atoms_list))+" does not exist. Please check your atoms indices")
        else:
            atoms_list = range(1, len(all_atoms_total_mulliken[0])+1)

        # Get angular contributions for projected DOS
        if args.angular:
            angular_list = args.angular
        else:
            angular_list = "a"

        for angular in angular_list:
            # Store this in case we are interpolating
            if len(total_mulliken) > 0:
                previous_total_mulliken = total_mulliken[:][:]
             
            # Catch erroneous values
            if angular not in angular_momenta:
                exit_error("Unknown option for angular momentum: "+angular+"\nPlease choose from spdf only")

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
                    #print electrons_on_this_specific_atom, normalising_constant, spin_counter

            # As above, calculate values and normalise. Could just copy from above but not so expensive.
            x, plot_values_mulliken, spin_counter = calculate_plot_values(min_point,max_point,args.polarised,energies,k_weights)
            normalise_plot_values(plot_values_mulliken, spin_counter, normalising_constant)

            # Add results to previous data
            if len(total_mulliken) > 0:
                total_mulliken = [ old + new for old, new in zip(total_mulliken, plot_values_mulliken) ]
            else:
                total_mulliken = plot_values_mulliken[:][:]
            
            if len(previous_total_mulliken) == 0:
                previous_total_mulliken = [[0] * len(tm) for tm in total_mulliken] 

            # Check the angular contribution is large enough to justify plotting
            total_change = 0.0
            for j in range(spin_counter):
                total_change = abs(sum(total_mulliken[j]) - sum(previous_total_mulliken[j]))

                # Arbitrary threshold of half an electron:
                if total_change > 0.5:
                    print("Plotting ", angular, " angular momenta in ", colors[angular_momenta[angular]-1], " for atoms ", atoms_list)
                    plt.plot(x,total_mulliken[j],lw ='2', color=colors[angular_momenta[angular]-1], label=angular, ls=line_types[k])
                    if args.interpolate:
                        plt.fill_between(x, previous_total_mulliken[j], total_mulliken[j],lw ='0', facecolor=colors[angular_momenta[angular]-1], label=file, interpolate=True)

    # Put this at the end so it covers everything else and shows the outline of the DOS correctly
    for j in range(spin_counter):
        plt.plot(x,total_plot_values[j],lw ='2', color='black', label=file, ls=line_types[k])

    # Work to rescale axes
    ymax = 1.1 
    ymin = 0
    if args.polarised:
        ymin = -ymax
    plt.ylim(ymin, ymax)
    plt.xlim(min_point+10, max_point-10)

    # HOMO
    plt.axvline(x=homo, ymin=-100, ymax=100, color='black', lw='2', ls=line_types[k]) # MFI

    # 0
    if args.polarised:
        plt.axhline(y=0, xmin=-100, xmax=100, color='black', lw='2')

    # Label, hacked
    if args.letters:
        plt.text(max_point-16, 0.95, labels[i], verticalalignment='top')

if args.align_to_homo:
    plt.xlabel(r'$\epsilon - \epsilon_{HOMO}$ (eV)')
else:
    plt.xlabel(r'$\epsilon$ (eV)')

# Display the graphs
print("done")
plt.show()
