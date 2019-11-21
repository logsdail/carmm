import sys
import pylab
import pickle
import os

pylab.rcParams.update({'font.size': 32})

input_files = sys.argv[1:]

if len(input_files) < 1:
    print "Error: Input file not found"
    sys.exit(1)

count = 0
width = None

for count in range(len(input_files)):
    filename = input_files[count]

    ef, energy, overall_dos = pickle.load(open(filename))
    #pylab.subplot(len(input_files), 1, count)

    if filename == "dos.spin1":
        overall_dos *= 2

    if os.path.exists(filename+".2"):
        pylab.plot(energy - ef[0], overall_dos)
        ef_2, energy_2, overall_dos_2 = pickle.load(open(filename+".2"))
        pylab.plot(energy_2 - ef_2[1], overall_dos_2)
        pylab.legend(('up', 'down'), loc='upper left')
    else:
        pylab.plot(energy - ef, overall_dos)
    
    pylab.xlabel(r'$\epsilon - \epsilon_F \ \rm{(eV)}$')
    pylab.ylabel('Density of States (1/eV)')

pylab.legend(input_files, loc='upper left')

pylab.axvline(x=0, ymin=0, ymax=1, ls='--', color='black', lw='1')

pylab.show()
