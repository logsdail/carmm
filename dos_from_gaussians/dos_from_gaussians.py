import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import sys
import pylab

pylab.rcParams.update({'font.size': 32})
input_files = sys.argv[1:]
energies = [[0 for x in xrange(0)] for x in xrange(len(input_files))]
colors = ['red','blue','violet','orange','blue','yellow','red','green','indigo']
interpolate_colors = ['Gold', 'MediumOrchid','LightCoral', 'LightGreen','LightPink','LightSkyBlue']

if len(input_files) < 1:
    print "Error: Input file not found"
    sys.exit(1)

ymax = 0

for i in range(len(input_files)):
    file = input_files[i]
    with open(file, 'r') as file_input:
        input_lines = file_input.readlines()
    for word in input_lines:
        energies[i].append(float(word))

    points = 1000
    values = [0]*points
    min_point = -18
    max_point = 5

    for mean in energies[i]:
        # Adjust specifically for PBC to match Ef
        #if i == 0:
                  #   VBM     + Vacuum Ef
            #mean += (-8.73891 + 8.4889732)
            #mean -= 0.43
        #    continue
        if mean > min_point and mean < max_point:
            variance = 0.05 
            sigma = np.sqrt(variance)
            x = np.linspace(min_point,max_point,points)
            values += mlab.normpdf(x,mean,sigma)

    # Normalise
    max_values = max(values)
    values /= max_values

    plt.plot(x,values,lw ='4', color=colors[i], label=input_files[i])

#    if i > 0:
#        plt.plot(x,values,lw ='4', color=colors[i], label=input_files[i])
#    else:
#        plt.plot(x,values,lw ='4', color='grey', label=input_files[i])

    # Work to rescale axes
    ycurrent = max(values)+0.1*max(values)
    if ycurrent > ymax:
        ymax = ycurrent
        pylab.axis([min_point, max_point, 0, ymax])

    # Interpolate beneath curves
    temp_values = [0]*points
    int_counter = 0
   
#    if i > 0:
#        continue 
    for j in range(len(values)):
    # Temp limit to highlight occupied states
        if int_counter < 2:
            if values[j] > 0.01:
                temp_values[j] = values[j]
                if j == len(values)-1:
                    pylab.fill_between(x, temp_values, y2=0, where=None, facecolor=interpolate_colors[int_counter], lw=0.0, interpolate=True, alpha=0.85)
            else:
                if j > 0:
                    if values[j-1] > 0.01:
                        pylab.fill_between(x, temp_values, y2=0, where=None, facecolor=interpolate_colors[int_counter], lw=0.0, interpolate=True, alpha=0.85)
                        temp_values = [0]*points
                        int_counter = int_counter + 1
#    else:
        # To deal with PBC
#    pylab.fill_between(x, values, y2=0, where=None, facecolor='white', lw=0.0, interpolate=True)
#    pylab.fill_between(x, values, y2=0, where=None, facecolor='grey', lw=0.0, interpolate=True, alpha=0.5)

pylab.rcParams.update({'mathtext.fontset': 'stix'})

plt.xlabel(r'$\epsilon$ (eV)')
plt.ylabel('Density of States (1/eV)')
# E_{f}
plt.axvline(x=-.84889732E+01, ymin=0, ymax=ymax, ls='--', color='black', lw='4') # MFI
#plt.axvline(x=-7.9, ymin=0, ymax=ymax, ls='--', color='black', lw='4')
#plt.legend(loc='upper right')

plt.show()

