#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
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
args =  parser.parse_args()

labels = ['(a)', '(b)', '(c)', '(d)', '(e)']
colors = ['red','blue','violet','orange','blue','yellow','red','green','indigo']
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

        energies = [0 for x in xrange(0)] 
        occupancy = [0 for x in xrange(0)] 
        homo = -999999.9
        lumo = 999999.9 

        file = input_files[i]
        with open(file, 'r') as file_input:
            input_lines = file_input.readlines()

        count=1

        #print file

        for sentence in input_lines:
            #print sentence
            split_sentence = sentence.split()
            #print split_sentence
            for word in split_sentence:
                if count%2 != 0:
                   #int(float(word))
                    energies.append(float(word))
                    #print "energy ", word
                else:
                    occupancy.append(float(word))
                    #print sentence
                count +=1
            #print count

        # Dirty hack if the file doesn't contain any results
        if len(energies) == 0:
            break 

        for j in range(len(occupancy)-1):
            #print occupancy[j]
            # Check if the occupancy of this state is more than 1, and less than one in the next state
            if occupancy[j] > 0.5 and occupancy[j+1] < 0.5:
                # Check if we already have a HOMO value that is higher, e.g. spin-polarised systems
                if energies[j] > homo:
                    homo = energies[j]
                if energies[j+1] < lumo:
                    lumo = energies[j+1]
                print "XC: ", xc,"; HOMO: ", homo, "; LUMO: ", lumo, "HOMO-LUMO: ", homo-lumo

        # Align everything to the homo?
        if args.align_to_homo:
            energies[:] = [y-homo for y in energies] 
            homo = 0.0
            plt.xlabel(r'$\epsilon - \epsilon_{HOMO}$ (eV)')

        points = 1000
        plot_values = [[0.0]*points]*2 #spin-up and spin-down
        min_point = homo - 20 
        max_point = homo + 20
        x = np.linspace(min_point,max_point,points)

        last_energy = energies[0]
	spin_counter = 1 

        #print xc, spin_counter, last_energy
        for energy in energies:
	    if args.split_updown and last_energy > energy:
	        spin_counter += 1
                #print xc, spin_counter, last_energy, energy
                  
	    last_energy = energy
            if energy > min_point and energy < max_point:
                variance = 0.02
                sigma = np.sqrt(variance)
		plot_values[spin_counter-1] += mlab.normpdf(x,energy,sigma)

	for j in range(spin_counter):
	    # Normalise
            max_value = max(plot_values[j])
            #plot_values /= max_value
            #No longer works so manually normalising over the array
            prefactor = 1
            if j > 0: 
	        prefactor = -1
            plot_values[j][:] = [prefactor*y/max_value for y in plot_values[j]]

        #Subplots
        ax1 = plt.subplot(len(input_files),1,1)
        if i > 0:
            ax1 = plt.subplot(len(input_files),1,1+i,sharex=ax1)
        # Hack to set the y label where I want it
        ax1.set_yticklabels([])
        if i == len(input_files)/2:
            plt.ylabel('Density of States (1/eV)')

	for j in range(spin_counter):
            #print j
            #Line plot
            if args.lines or not args.interpolate:
                plt.plot(x,plot_values[j],lw ='4', color=colors[i], label=input_files[i], ls=line_types[k])
            else:
                plt.plot(x,plot_values[j],lw ='2', color='black', label=input_files[i], ls=line_types[k])
            #Interpolate plot
            if args.interpolate:
                plt.fill_between(x, 0, plot_values[j],lw ='0', facecolor=colors[i], label=input_files[i], interpolate=True, alpha=0.5)

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
        plt.axvline(x=homo, ymin=-100, ymax=100, color=colors[i], lw='2', ls=line_types[k]) # MFI

        # 0
        if args.split_updown:
            plt.axhline(y=0, xmin=-100, xmax=100, color='black', lw='2')

        # Label, hacked
        if args.letters:
            plt.text(-4, 0.95, labels[i], verticalalignment='top')

#plt.legend(loc='upper right')

# Make sure this label is on the bottom graph
#ax1 = plt.subplot(len(input_files),1,len(input_files),sharex=ax1)
# Hack to set the y label where I want it
#ax1.set_yticklabels([])
plt.xlabel(r'$\epsilon$ (eV)')
#plt.ylabel('Density of States (1/eV)')
plt.show()

