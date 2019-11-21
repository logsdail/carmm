import sys
import pylab
import os
import pickle
import numpy
from ase.data import chemical_symbols

pylab.rcParams.update({'font.size': 32})
colours = ['red','blue','green','orange','violet','yellow','indigo']
#colours = ['red','orange','grey'] 

input_files = sys.argv[1:]
if len(input_files) < 1:
    print "Error: Input file not found"
    sys.exit(1)

ef = [[] for x in xrange(len(input_files))]
elements = [[] for x in xrange(len(input_files))]
energy = [[] for x in xrange(len(input_files))]
dos = [[] for x in xrange(len(input_files))]
dos_occ = [[] for x in xrange(len(input_files))]
dos_occ_com = [[] for x in xrange(len(input_files))]
dos_unocc = [[] for x in xrange(len(input_files))]
dos_unocc_com = [[] for x in xrange(len(input_files))]

for ifile in range(len(input_files)):
    filename = input_files[ifile]

    # Padding
    print " ======================= ", filename, " ======================= "
    print ""

    factor = 2
    if os.path.exists(filename+".spin2") or os.path.exists(filename+".spin1.2"):
        factor = 1

#    for element in chemical_symbols:
    for element in ['core','face','shell','edge']:

        if os.path.exists(filename+"."+element):
            temp_ef, temp_energy, temp_dos = pickle.load(open(filename+"."+element))
            print element, " found for: ", filename, "; energy scale from ", temp_energy[1], " to ", temp_energy[-1]

            # Scale up in the case of spin-paired calculations
            temp_dos *= factor
            
            # Append information to various arrays
            ef[ifile].append(temp_ef)
            elements[ifile].append(element)
            energy[ifile].append(temp_energy)
            dos[ifile].append(temp_dos)

            dos_occ[ifile].append([0]*len(temp_dos))
            dos_occ_com[ifile].append(0)
            dos_unocc[ifile].append([0]*len(temp_dos))
            dos_unocc_com[ifile].append(0)

	    # Sort out occupied and unoccupied dos
            for ienergy in range(len(temp_energy)):
                if temp_ef > temp_energy[ienergy]:
                    dos_occ_com[ifile][len(dos_occ_com[ifile])-1] += temp_energy[ienergy] * temp_dos[ienergy]
                    dos_occ[ifile][len(dos_unocc[ifile])-1][ienergy] = temp_dos[ienergy] 
                else:
                    dos_unocc_com[ifile][len(dos_unocc_com[ifile])-1] += temp_energy[ienergy]*temp_dos[ienergy]
                    dos_unocc[ifile][len(dos_unocc[ifile])-1][ienergy] = temp_dos[ienergy]

            dos_occ_com[ifile][len(dos_occ_com[ifile])-1]  /= sum(dos_occ[ifile][len(dos_occ[ifile])-1]) 
            dos_unocc_com[ifile][len(dos_unocc_com[ifile])-1] /= sum(dos_unocc[ifile][len(dos_unocc[ifile])-1])
        
            print "Centre of Mass: Occupied: ", dos_occ_com[ifile][len(dos_occ_com[ifile])-1], "; Unoccupied: ", dos_unocc_com[ifile][len(dos_unocc_com[ifile])-1]      
            print "" 
            print "Trapezium rule: Occupied States + Unoccupied States = Total States" 
            print numpy.trapz(dos_occ[ifile][len(dos_occ[ifile])-1], temp_energy), " + ", numpy.trapz(dos_unocc[ifile][len(dos_unocc[ifile])-1], temp_energy), " = ", numpy.trapz(temp_dos, temp_energy)
            print ""

    pylab.subplot(len(input_files), 1, ifile)

# Interpolated filled plots STRONG COLOURS
    colour_counter = 0

    if len(ef[ifile]) > 1:
        # We have more than one speces. Calculate a combined COM for the DOS
        total_dos_occ = dos_occ[ifile][0]
        total_dos_unocc = dos_unocc[ifile][0]
        total_dos_occ_com = dos_occ_com[ifile][0]*sum(dos_occ[ifile][0])
        total_dos_unocc_com = dos_unocc_com[ifile][0]*sum(dos_unocc[ifile][0])

        for ielement in range(1,len(ef[ifile])):
            total_dos_occ += dos_occ[ifile][ielement]
            total_dos_unocc += dos_unocc[ifile][ielement]
            total_dos_occ_com += dos_occ_com[ifile][ielement]*sum(dos_occ[ifile][ielement])
            total_dos_unocc_com += dos_unocc_com[ifile][ielement]*sum(dos_unocc[ifile][ielement])

        print "Centre of Mass Total DOS: Occupied: ", total_dos_occ_com/sum(total_dos_occ), "; Unoccupied: ", total_dos_unocc_com/sum(total_dos_unocc)
        print ""

        # Now to superimpose these dos profiles using their individual scales
        # Energy scales are equivlent (fortunately!)

        previous = [0]*len(dos[ifile][0])
        current = [0]*len(dos[ifile][0])
        print "Ef for this series (check they are consistent!)"

        for ielement in range(len(ef[ifile])):
            current += dos[ifile][ielement]
            print elements[ifile][ielement], ef[ifile][ielement], colours[colour_counter]
            pylab.fill_between(energy[ifile][0], current, previous, facecolor=colours[colour_counter], lw=0.0, interpolate=True, label=elements[ifile][ielement])
            # pylab.plot(energy[ifile][0], current, color=colours[colour_counter], lw=2.0, label=elements[ifile][ielement]) 
            colour_counter += 1
            previous = numpy.copy(current)

        print ""       

    else:
        pylab.fill_between(energy[ifile][0], dos[ifile][0], y2=0, where=None, facecolor=colours[colour_counter], lw=0.0, interpolate=True, label=elements[ifile][0])
#        pylab.fill_between(energy[ifile][0], dos[ifile][0], y2=0, where=None, facecolor='orange', lw=0.0, interpolate=True, label=elements[ifile][0])

    pylab.legend()

    xmin = -15
    xmax = 0 
    ymin = 0
    ymax = 500
    pylab.axis([xmin, xmax, ymin, ymax]) 
    pylab.yticks( range(0) )
    pylab.text((xmin+xmax)/1.1,(ymax-ymin)/1.4,filename)
    pylab.axvline(x=ef[ifile][0], ymin=0, ymax=1, ls='--', color='black', lw='1.0')

    if ifile == 0:
        pylab.xlabel(r'$\epsilon$')
        pylab.xticks(range(xmin,xmax,5))
    else:
        pylab.xticks( range(0) )

    if ifile == ((len(input_files)/2)):
        pylab.ylabel('Density of States (1/eV)')
pylab.show()
