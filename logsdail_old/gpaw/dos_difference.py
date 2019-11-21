import sys
import pylab
from gpaw import GPAW
import matplotlib
import pickle
#matplotlib.use('Agg')
import os

pylab.rcParams.update({'font.size': 32})

files_length = 2 
files = ['']*files_length

for i in range(files_length):
    files[i] = sys.argv[i+1]

count = 0
width = None

#pylab.subplot(2, 1, 1)

for filename in files:
    count += 1

    ef, energy, overall_dos = pickle.load(open(filename))

    print energy[1], energy[-1]

    #pylab.subplot(len(files), 1, count)
    #pylab.plot(energy - ef, overall_dos)

    if os.path.exists(filename+".2"):
        ef_2, energy_2, overall_dos_2 = pickle.load(open(filename+".2"))
        overall_dos += overall_dos_2
       #pylab.legend(('up', 'down'), loc='upper left')
    else:
        overall_dos *= 2

    if count == 1:
        pylab.plot(energy, overall_dos, lw ='2', color='Green') 
        overall_dos1 = overall_dos
        ef1 = ef
    else:
        pylab.plot(energy, overall_dos, lw ='2', ls='--', color='Blue')
        overall_dos2 = overall_dos
        ef2 = ef

pylab.xlabel(r'$\epsilon$', fontsize=32)
pylab.ylabel('Density of States (1/eV)',fontsize=32)

#ax = pylab.axes()
#ax.plot(pylab.arange(10))
#xlabels = ax.get_xticklabels()
#xlabel0 = xlabels[0]                              # one of the xtick labels
#xlabel0.get_fontsize()
#xlabel0.set_fontsize(20) 

overall_dos = overall_dos1 - overall_dos2

#pylab.subplot(2, 1, 2)
pylab.axvline(x=ef1,ymin=0, ymax=1, ls='--', color='black', lw='2')
pylab.axvline(x=ef2, ymin=0, ymax=1, ls='--', color='grey', lw='2')
#pylab.fill_between(energy, overall_dos, y2=0, where=None, facecolor='Red', lw=0, interpolate=True)
pylab.axhline(y=0, xmin=-100, xmax=100, color='black', lw='2')

xmin = -15
xmax = 0
ymin = -200
ymax = 500 
pylab.axis([xmin, xmax, ymin, ymax])  
#plt.yticks( range(0) )

pylab.show()
