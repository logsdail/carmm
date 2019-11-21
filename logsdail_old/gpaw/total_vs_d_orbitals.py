import sys
import pylab
import os
import pickle
import numpy

pylab.rcParams.update({'font.size': 32})

files = sys.argv[1:]
if len(files) < 1:
    print "Error: Input file not found"
    sys.exit(1)

files_length = len(files)

count = 0
width = 0.2
labels = ['(a)','(b)','(c)','(d)','(e)'] 
labels_2 = ['Au','Ag','Au@Ag','Ag@Au','Pd@Au (Ico)']
labels_3 = ['N = 55', 'N = 147', 'N = 309']
labels_4 = ['Ico','Dec','Cub']

for filename in files:
    count += 1

# Calculate overall pdos

    overall_dos = []
    factor = 2 
#
    ef, energy, overall_dos = pickle.load(open(filename+".spin1")) 
    if os.path.exists(filename+".spin2") or os.path.exists(filename+".spin1.2"):
        if os.path.exists(filename+".spin2"):
            ef_2, energy_2, overall_dos_2 = pickle.load(open(filename+".spin2"))
        else:
            ef_2, energy_2, overall_dos_2 = pickle.load(open(filename+".spin1.2"))

        overall_dos += overall_dos_2
        factor = 1 
    else:
        overall_dos *= factor

    # Padding
    print " ======================= ", filename, " ======================= "
    print ""
    print "Trapezium rule: Total DOS States = ", numpy.trapz(overall_dos, energy)
    print ""

    if os.path.exists(filename+".core.d"):
        ef, au_energy, au_dos = pickle.load(open(filename+".core.d"))
        print "Core found for: ", filename, "; energy scale from ", au_energy[1], " to ", au_energy[-1]

        # Scale up in the case of spin-paired calculations
        au_dos *= factor

        au_dos_occ = [0]*len(au_dos)
        au_dos_occ_com = 0
        au_dos_occ_sum = 0
        au_dos_unocc = [0]*len(au_dos)
        au_dos_unocc_com = 0
        au_dos_unocc_sum = 0
        for i in range(len(au_energy)):
            if ef > au_energy[i]:
                au_dos_occ[i] = au_dos[i]
                au_dos_occ_sum += au_dos[i]
                au_dos_occ_com += au_dos[i] * au_energy[i]
            else:
                au_dos_unocc[i] = au_dos[i]
                au_dos_unocc_sum += au_dos[i]
                au_dos_unocc_com += au_dos[i] * au_energy[i]
        
        print "Centre of Mass: Occupied: ", au_dos_occ_com/au_dos_occ_sum, "; Unoccupied: ", au_dos_unocc_com/au_dos_unocc_sum       
        print "" 
        print "Trapezium rule: Occupied States + Unoccupied States = Total States" 
        print numpy.trapz(au_dos_occ, au_energy), " + ", numpy.trapz(au_dos_unocc, au_energy), " = ", numpy.trapz(au_dos, au_energy)
        print ""

    if os.path.exists(filename+".shell.d"):
        ef, pd_energy, pd_dos = pickle.load(open(filename+".shell.d"))
        print "Shell found for: ", filename, "; energy scale from ", pd_energy[1], " to ", pd_energy[-1]

        # Scale up in the case of spin-paired calculations
        pd_dos *= factor

        pd_dos_occ = [0]*len(pd_dos)
        pd_dos_occ_com = 0
        pd_dos_occ_sum = 0
        pd_dos_unocc = [0]*len(pd_dos)
        pd_dos_unocc_com = 0
        pd_dos_unocc_sum = 0
        for i in range(len(pd_energy)):
            if ef > pd_energy[i]:
                pd_dos_occ[i] = pd_dos[i]
                pd_dos_occ_sum += pd_dos[i]
                pd_dos_occ_com += pd_dos[i] * pd_energy[i]
            else:
                pd_dos_unocc[i] = pd_dos[i]
                pd_dos_unocc_sum += pd_dos[i]
                pd_dos_unocc_com += pd_dos[i] * pd_energy[i]

        print "Centre of Mass: Occupied: ", pd_dos_occ_com/pd_dos_occ_sum, "; Unoccupied: ", pd_dos_unocc_com/pd_dos_unocc_sum
        print ""
        print "Trapezium rule: Occupied States + Unoccupied States = Total States"
        print numpy.trapz(pd_dos_occ, pd_energy), " + ", numpy.trapz(pd_dos_unocc, pd_energy), " = ", numpy.trapz(pd_dos, pd_energy)
        print ""

    pylab.subplot(len(files), 1, count)

# Interpolated filled plots STRONG COLOURS
    pylab.plot(energy, overall_dos, color='grey', lw='2.0')

    print "Fermi energy          : ", ef

    array_index = numpy.where(energy==(min(energy, key=lambda x:abs(x-ef))))[0][0]

    if energy[array_index] < ef:
        array_index=array_index+1

    print "Closest to it         : ", energy[array_index], " at index ", array_index
    # print au_energy[array_index], pd_energy[array_index]
    print "No. of states         : ", overall_dos[array_index]     

    if os.path.exists(filename+".core.d") and os.path.exists(filename+".shell.d"):

        print "No. of core d-states  : ", au_dos[array_index]
        print "No. of shell d-states : ", pd_dos[array_index]

        total_dos = au_dos + pd_dos
        total_dos_occ_com = au_dos_occ_com + pd_dos_occ_com 
        total_dos_occ_sum = au_dos_occ_sum + pd_dos_occ_sum
        total_dos_unocc_com = au_dos_unocc_com + pd_dos_unocc_com
        total_dos_unocc_sum = au_dos_unocc_sum + pd_dos_unocc_sum

        print "Centre of Mass Total DOS: Occupied: ", total_dos_occ_com/total_dos_occ_sum, "; Unoccupied: ", total_dos_unocc_com/total_dos_unocc_sum
        print ""

        pylab.plot(au_energy, au_dos, color='red', lw='2.0') #core
        pylab.plot(pd_energy, pd_dos, color='blue', lw='2.0') #core
        #pylab.fill_between(au_energy, total_dos, au_dos, where=None, facecolor='orange', lw=0.0, interpolate=True)
        #pylab.fill_between(pd_energy, au_dos, y2=0, where=None, facecolor='blue', lw=0.0, interpolate=True)

#    else:
#        pylab.fill_between(energy - ef, overall_dos, y2=0, where=None, facecolor='grey', lw=0.0, interpolate=True)

#pylab.legend()
    if count == files_length:
        pylab.xlabel(r'$\epsilon$')
        pylab.xticks(range(-12,2,2))
    else:
        pylab.xticks( range(0) )

    xmin = -12
    xmax = -2 
    ymin = 0
    ymax = 700
    pylab.axis([xmin, xmax, ymin, ymax]) 
    pylab.yticks( range(0) )
#    pylab.text(-10,500,filename)
    pylab.text(xmin+1,ymax-200,labels_2[count-1])
    pylab.axvline(x=ef, ymin=0, ymax=1, ls='--', color='black', lw='2.0')

    if count == ((files_length/2)+1):
        pylab.ylabel('Density of States (1/eV)')
pylab.show()
