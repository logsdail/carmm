import sys
import pylab
from gpaw import GPAW
import matplotlib
import pickle
#matplotlib.use('Agg')
import os
import numpy

pylab.rcParams.update({'font.size': 32})

files_length = len(sys.argv)-1 
files = ['']*files_length
color = ['Red','Green','Blue','Orange']
linetypes = ['-','--','.',':']

for i in range(files_length):
    files[i] = sys.argv[i+1]

count = 0
width = None
max_dos = 0.0
min_energy = -99999.0
max_energy = 99999.0

#pylab.subplot(2, 1, 1)

for filename in files:

    ef, energy, overall_dos = pickle.load(open(filename))
    if count < 2:
        overall_dos /= 2.0
    if count > 2:
        overall_dos /= 2.0

    #pylab.subplot(len(files), 1, count)
    #pylab.plot(energy - ef, overall_dos)

    if os.path.exists(filename[:-1]+"2"):
        ef_2, energy_2, overall_dos_2 = pickle.load(open(filename[:-1]+"2"))
        overall_dos += overall_dos_2
        #pylab.legend(('up', 'down'), loc='upper left')

    overall_dos += count*100

    pylab.plot(energy, overall_dos, lw ='2', color=color[count]) #ls=linetypes[count] 

    max_dos = max(max_dos,numpy.amax(overall_dos))
    min_energy = max(min_energy,min(energy))
    max_energy = min(max_energy,max(energy))

    print ef 
    pylab.axvline(x=ef, ymin=0, ymax=1, ls='--', color=color[count], lw='2')
    count += 1

pylab.xlabel(r'$\epsilon \ \rm{(eV)}$', fontsize=32)
pylab.ylabel('Density of States (1/eV)',fontsize=32)

#ax = pylab.axes()
#ax.plot(pylab.arange(10))
#xlabels = ax.get_xticklabels()
#xlabel0 = xlabels[0]                              # one of the xtick labels
#xlabel0.get_fontsize()
#xlabel0.set_fontsize(20) 

#ef, energy, overall_dos1 = pickle.load(open(files[0]))
#ef, energy, overall_dos2 = pickle.load(open(files[1]))
#overall_dos = overall_dos1 - overall_dos2

#pylab.subplot(2, 1, 2)
#pylab.axvline(x=0, ymin=0, ymax=1, ls='--', color='black', lw='2')
#pylab.fill_between(energy - ef, overall_dos, y2=0, where=None, facecolor='Red', lw=0, interpolate=True)
#pylab.axhline(y=0, xmin=-100, xmax=100, color='black', lw='2')

xmin = min_energy + 1
xmax = max_energy - 1
ymin = 0
ymax = max_dos + 10
pylab.axis([xmin, xmax, ymin, ymax])  
pylab.yticks( range(0) )

pylab.show()
