import sys
import pickle
import numpy as np
import pylab as plt
#from gpaw import GPAW
#import matplotlib
#matplotlib.use('Agg')
import os

#from gpaw.utilities.dos import print_projectors
#print_projectors('Au')

plt.rcParams.update({'font.size': 32})

input_files = sys.argv[1:]
files_length = len(input_files)

if len(input_files) < 1:
    print "Error: Input file not found"
    sys.exit(1)

#### COUNTERS ####
total_sum = [0]*files_length
total_com = [0]*files_length
s_dos_occ_com = [0]*files_length
s_dos_occ_sum = [0]*files_length
p_dos_occ_com = [0]*files_length
p_dos_occ_sum = [0]*files_length
d_dos_occ_com = [0]*files_length
d_dos_occ_sum = [0]*files_length
f_dos_occ_com = [0]*files_length
f_dos_occ_sum = [0]*files_length

for ifile in range(len(input_files)):
    input = input_files[ifile]

    # Padding
    print " ======================= ", input, " ======================= "
    print ""


    #ef, e_spdf, p_spdf = pickle.load(open(input))

    # Total range is ~5 to -20.
    # And we want a resolution of 201pts for smaller range
    # Therefore increase to deal with this

    s = False
    p = False
    d = False
    f = False
    npts = 0

    if os.path.exists(input+".s"):
        s = True
        ef, e_s, s_dos = pickle.load(open(input+".s"))
        npts = len(e_s)
    if os.path.exists(input+".p"):   
        p = True
        ef, e_p, p_dos = pickle.load(open(input+".p"))
    if os.path.exists(input+".d"):   
        d = True
        ef, e_d, d_dos = pickle.load(open(input+".d"))
    if os.path.exists(input+".f"):   
        f = True 
        ef, e_f, f_dos = pickle.load(open(input+".f"))

    if ifile < 2:
        s_dos /= 2
        p_dos /= 2
        d_dos /= 2
        f_dos /= 2 

    ef_point = 0

    for i in range(npts):
        if ef <= e_s[i] and ef_point == 0:
            ef_point = i
            #print e_s[i], ef, ef_point
    dos = [0]*ef_point
    e = [0]*ef_point

    if s:
    ## S

        s_dos_occ = [0]*len(s_dos)
#        s_dos_occ_com = 0
#        s_dos_occ_sum = 0
        s_dos_unocc = [0]*len(s_dos)
        s_dos_unocc_com = 0
        s_dos_unocc_sum = 0

        for i in range(npts):
            if ef > e_s[i]:
                s_dos_occ[i] = s_dos[i]
                s_dos_occ_sum[ifile] += s_dos[i]
                s_dos_occ_com[ifile] += s_dos[i] * e_s[i]
            else:
                s_dos_unocc[i] = s_dos[i]
                s_dos_unocc_sum += s_dos[i]
                s_dos_unocc_com += s_dos[i] * e_s[i]

        for i in range(ef_point):
            e[i] = e_s[i]
            dos[i] = s_dos[i] 

        print 's:   occupied: sum:', s_dos_occ_sum[ifile], 'average:', s_dos_occ_com[ifile]/s_dos_occ_sum[ifile], 'integral:', np.trapz(s_dos_occ, e_s)
#        print 's: unoccupied: sum:', s_dos_unocc_sum, 'average:', s_dos_unocc_com/s_dos_unocc_sum, 'integral:', np.trapz(s_dos_unocc, e_s)
    
    if p:
    ## P

        p_dos_occ = [0]*len(p_dos)
#        p_dos_occ_com = 0
#        p_dos_occ_sum = 0
        p_dos_unocc = [0]*len(p_dos)
        p_dos_unocc_com = 0
        p_dos_unocc_sum = 0

        for i in range(npts):
            if ef > e_p[i]:
                p_dos_occ[i] = p_dos[i]
                p_dos_occ_sum[ifile] += p_dos[i]
                p_dos_occ_com[ifile] += p_dos[i] * e_p[i]
            else:
                p_dos_unocc[i] = p_dos[i]
                p_dos_unocc_sum += p_dos[i]
                p_dos_unocc_com += p_dos[i] * e_p[i]

        print 'p:   occupied: sum:', p_dos_occ_sum[ifile], 'average:', p_dos_occ_com[ifile]/p_dos_occ_sum[ifile], 'integral:', np.trapz(p_dos_occ, e_p)
#        print 'p: unoccupied: sum:', p_dos_unocc_sum, 'average:', p_dos_unocc_com/p_dos_unocc_sum, 'integral:', np.trapz(p_dos_unocc, e_p)    

        for i in range(ef_point):
            e[i] = e_p[i]
            dos[i] = p_dos[i]

    ## ADD
 
    p_dos_filled = s_dos + p_dos

    if d:
    ## D

        d_dos_occ = [0]*len(d_dos)
#        d_dos_occ_com = 0
#        d_dos_occ_sum = 0
        d_dos_unocc = [0]*len(d_dos)
        d_dos_unocc_com = 0
        d_dos_unocc_sum = 0

        for i in range(npts):
            if ef > e_d[i]:
                d_dos_occ[i] = d_dos[i]
                d_dos_occ_sum[ifile] += d_dos[i]
                d_dos_occ_com[ifile] += d_dos[i] * e_d[i]
            else:
                d_dos_unocc[i] = d_dos[i]
                d_dos_unocc_sum += d_dos[i]
                d_dos_unocc_com += d_dos[i] * e_d[i]

        print 'd:   occupied: sum:', d_dos_occ_sum[ifile], 'average:', d_dos_occ_com[ifile]/d_dos_occ_sum[ifile], 'integral:', np.trapz(d_dos_occ, e_d)
#        print 'd: unoccupied: sum:', d_dos_unocc_sum, 'average:', d_dos_unocc_com/d_dos_unocc_sum, 'integral:', np.trapz(d_dos_unocc, e_d)

        for i in range(ef_point):
            e[i] = e_d[i]
            dos[i] = d_dos[i]

    ## ADD

    d_dos_filled = p_dos_filled + d_dos

    if f:
    ## F

        f_dos_occ = [0]*len(f_dos)
#        f_dos_occ_com = 0
#        f_dos_occ_sum = 0
        f_dos_unocc = [0]*len(f_dos)
        f_dos_unocc_com = 0
        f_dos_unocc_sum = 0

        for i in range(npts):
            if ef > e_f[i]:
                f_dos_occ[i] = f_dos[i]
                f_dos_occ_sum[ifile] += f_dos[i]
                f_dos_occ_com[ifile] += f_dos[i] * e_f[i]
            else:
                f_dos_unocc[i] = f_dos[i]
                f_dos_unocc_sum += f_dos[i]
                f_dos_unocc_com += f_dos[i] * e_f[i]
        
        if f_dos_occ_sum[ifile] > 0:
            print 'f:   occupied: sum:', f_dos_occ_sum[ifile], 'average:', f_dos_occ_com[ifile]/f_dos_occ_sum[ifile], 'integral:', np.trapz(f_dos_occ, e_f)
#        print 'f: unoccupied: sum:', f_dos_unocc_sum, 'average:', f_dos_unocc_com/f_dos_unocc_sum, 'integral:', np.trapz(f_dos_unocc, e_f)

        for i in range(ef_point):
            e[i] = e_f[i]
            dos[i] = f_dos[i]

    ## ADD

    f_dos_filled = d_dos_filled + f_dos

    #plt.autumn()
    plt.subplot(files_length, 1, ifile)

    total_dos_filled = []

    if s:
    #plt.plot(energy - ef, s_dos, color='green')
    #plt.fill_between(energy - ef, s_dos, y2=0, where=None, facecolor='green', edgecolor='green', interpolate=True)
    #plt.fill_between(energy - ef, s_dos, y2=0, where=None, facecolor='DarkRed', edgecolor='DarkRed', interpolate=True)
    #plt.fill_between(energy - ef, s_dos, y2=0, where=None, facecolor='LightBlue', lw=0, interpolate=True)
        plt.fill_between(e_s, s_dos, y2=0, where=None, facecolor='Blue', lw=0, interpolate=True)
    #plt.fill_between(energy - ef, s_dos, y2=0, where=None, interpolate=True)
        total_dos_filled = s_dos

    if p:
    #plt.plot(energy - ef, p_dos)
    #plt.plot(energy - ef, p_dos_filled, color='blue')
    #plt.fill_between(energy - ef, p_dos, y2=0, where=None, facecolor='blue', edgecolor='blue', interpolate=True)
    #plt.fill_between(energy - ef, p_dos_filled, s_dos, where=None, facecolor='red', edgecolor='red', interpolate=True)
    #plt.fill_between(energy - ef, p_dos_filled, s_dos, where=None, facecolor='DarkGreen', lw=0, interpolate=True)
        plt.fill_between(e_p, p_dos_filled, s_dos, where=None, facecolor='Yellow', lw=0, interpolate=True)
    #plt.fill_between(energy - ef, p_dos_filled, s_dos, where=None, interpolate=True)
        total_dos_filled = p_dos_filled

    if d:
    #plt.plot(energy - ef, d_dos)
    #plt.plot(energy - ef, d_dos_filled, color='yellow')
    #plt.fill_between(energy - ef, d_dos, y2=0, where=None, facecolor='yellow', edgecolor='yellow', interpolate=True)
    #plt.fill_between(energy - ef, d_dos_filled, p_dos_filled, where=None, facecolor='orange', edgecolor='orange', interpolate=True)
    #plt.fill_between(energy - ef, d_dos_filled, p_dos_filled, where=None, facecolor='LightGreen', lw=0, interpolate=True)
        plt.fill_between(e_d, d_dos_filled, p_dos_filled, where=None, facecolor='Red', lw=0, interpolate=True)
    #plt.fill_between(energy - ef, d_dos_filled, p_dos_filled, where=None, interpolate=True)
        total_dos_filled = d_dos_filled

    if f:
    #plt.plot(energy - ef, f_dos)
    #plt.plot(energy - ef, f_dos_filled, color='red')
    #plt.fill_between(energy - ef, f_dos, y2=0, where=None, facecolor='red', edgecolor='red', interpolate=True)
        plt.fill_between(e_f, f_dos_filled, d_dos_filled, where=None, facecolor='Orange', lw=0, interpolate=True)
    #plt.fill_between(energy - ef, f_dos_filled, d_dos_filled, where=None, facecolor='LightGreen', edgecolor='LightGreen', interpolate=True)
    #plt.fill_between(energy - ef, f_dos_filled, d_dos_filled, where=None, interpolate=True)
        total_dos_filled = f_dos_filled

    print ""
    print "Fermi energy          : ", ef
    print ""

    array_index = np.where(e_s==(min(e_s, key=lambda x:abs(x-ef))))[0][0]

    if e_s[array_index] < ef:
        array_index=array_index+1

    print "Closest to it         : ", e_s[array_index], " at index ", array_index
    print "Total No. of states   : ", total_dos_filled[array_index]
    print "No. of d-states       : ", d_dos[array_index]
    print ""

#    plt.legend()
#    if count == files_length:
#    plt.xlabel(r'$\epsilon - \epsilon_F \ \rm{(eV)}$')
#    plt.ylabel('Density of States ')
    if ifile == 0:
        plt.xlabel(r'$\epsilon - \epsilon_F \ \rm{(eV)}$')
#        pylab.ylabel('Density of States (1/eV)')
    else:
        plt.xticks( range(0) )

    if ifile == (files_length - 1):
        plt.ylabel('Density of States (1/eV)')

    plt.tick_params(left=False,right=False,top=False)

#    xmin = min(e_s)
#    xmax = max(e_s)
    xmin = -12
    xmax = -2 
    ymin = 0
    ymax = 350 
    plt.axis([xmin, xmax, ymin, ymax])
    plt.yticks( range(0) )
    plt.text(-14,200,input,fontsize=12)
 
    plt.axvline(x=ef, ymin=0, ymax=1, ls='--', color='black', lw='2')
    
    plt.text((d_dos_occ_com[ifile]/d_dos_occ_sum[ifile])-0.05,0.55*ymax,r'$\bar{\epsilon}^{d}$')
    plt.axvline(x=d_dos_occ_com[ifile]/d_dos_occ_sum[ifile], ymin=0, ymax=0.5, ls='-', color='black', lw='2')

    total_sum[ifile] = s_dos_occ_sum[ifile]+d_dos_occ_sum[ifile]+p_dos_occ_sum[ifile]+f_dos_occ_sum[ifile]
    total_com[ifile] = s_dos_occ_com[ifile]+d_dos_occ_com[ifile]+p_dos_occ_com[ifile]+f_dos_occ_com[ifile]

    plt.text((total_com[ifile]/total_sum[ifile])-0.05,0.55*ymax,r'$\bar{\epsilon}$')
    plt.axvline(x=total_com[ifile]/total_sum[ifile], ymin=0, ymax=0.5, ls=':', color='black', lw='2')

    print 'tot: occupied: sum:', total_sum[ifile], 'average:', total_com[ifile]/total_sum[ifile]
#    print 'total: unoccupied: sum:', s_dos_unocc_sum+d_dos_unocc_sum+p_dos_unocc_sum+f_dos_unocc_sum, 'average:', (s_dos_unocc_com+d_dos_unocc_com+p_dos_unocc_com+f_dos_unocc_com)/(s_dos_unocc_sum+d_dos_unocc_sum+p_dos_unocc_sum+f_dos_unocc_sum)

total_s_occ_com = 0
total_s_occ_sum = 0
total_p_occ_com = 0
total_p_occ_sum = 0
total_d_occ_com = 0
total_d_occ_sum = 0
total_f_occ_com = 0
total_f_occ_sum = 0

for i in range(files_length):
    total_s_occ_com += s_dos_occ_com[i]
    total_s_occ_sum += s_dos_occ_sum[i]
    total_p_occ_com += p_dos_occ_com[i]
    total_p_occ_sum += p_dos_occ_sum[i]
    total_d_occ_com += d_dos_occ_com[i]
    total_d_occ_sum += d_dos_occ_sum[i]
    total_f_occ_com += f_dos_occ_com[i]
    total_f_occ_sum += f_dos_occ_sum[i]


print 'mean:'
print 's:   occupied: sum:', total_s_occ_sum, 'average:', total_s_occ_com/total_s_occ_sum
print 'p:   occupied: sum:', total_p_occ_sum, 'average:', total_p_occ_com/total_p_occ_sum
print 'd:   occupied: sum:', total_d_occ_sum, 'average:', total_d_occ_com/total_d_occ_sum
if total_f_occ_sum > 0:
    print 'f:   occupied: sum:', total_f_occ_sum, 'average:', total_f_occ_com/total_f_occ_sum
print 'tot: occupued: sum:', total_s_occ_sum+total_p_occ_sum+total_d_occ_sum+total_f_occ_sum, 'average:', (total_s_occ_com+total_p_occ_com+total_d_occ_com+total_f_occ_com)/(total_s_occ_sum+total_p_occ_sum+total_d_occ_sum+total_f_occ_sum)

if files_length == 2:
    print '',s_dos_occ_com[0]/s_dos_occ_sum[0], p_dos_occ_com[0]/p_dos_occ_sum[0], d_dos_occ_com[0]/d_dos_occ_sum[0], '-', total_com[0]/total_sum[0], s_dos_occ_com[1]/s_dos_occ_sum[1], p_dos_occ_com[1]/p_dos_occ_sum[1], d_dos_occ_com[1]/d_dos_occ_sum[1], '-', total_com[1]/total_sum[1],total_s_occ_com/total_s_occ_sum, total_p_occ_com/total_p_occ_sum, total_d_occ_com/total_d_occ_sum, '-', (total_s_occ_com+total_p_occ_com+total_d_occ_com+total_f_occ_com)/(total_s_occ_sum+total_p_occ_sum+total_d_occ_sum+total_f_occ_sum)

plt.show()
