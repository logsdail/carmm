def distance_distribution_function(model, bin_sampling):
    '''Returns a plot of the distribution of the distance between all atoms
    plot is currently a frequency vs distance. Current usage is for periodic solids
   TODO:
         - Only calculates RDF with respect the first atom - needs to be generalised:
            - what about for any other atom of interest?
            - what about averaging over all atoms, as per standard EXAFS?
              OB: amended distance_distribution to average over all atoms
         - Need to account for density

        '''

    from matplotlib import pyplot as plt
    from math import ceil
    import numpy as np
    import pylab as pl

    # get all distances in the model
    distances = model.get_all_distances(mic=True, vector=False)

    individual_lengths = []
    for i in range(len(distances)):
        for j in range(i+1, len(distances[i])):
            individual_lengths.append(distances[i][j])

    # plot these values as a histogram then as a line

    y, binEdges = np.histogram(individual_lengths, bins=ceil(max(individual_lengths) / bin_sampling))
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-')
    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('g(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Distribution Function', fontsize=15)
    pl.show()
    plt.show()



def radial_distribution_function(model, radius, position):
    '''Returns a plot of the distribution of the distance between each atom from atom_0.
    plot is currently a frequency vs distance. Current usage is for periodic solids
    This script will create a radius around a given model and calculate the distance of this new model
    TODO:  - Gaussian over the histogram
            AJL: What's the difference between this and difference_distribution_function?
    '''
    from matplotlib import pyplot as plt
    from carmm.build.cutout import cutout_sphere
    import numpy as np
    import pylab as pl
    #Read file or Atoms object

    if isinstance(model, str) is True:
        from ase.io import read
        model = read(model)
    # Create a variable which represents the amount of atoms in the system
    positions = model.get_positions()
    number_atoms_in_model = len(positions)
    # create a super cell of the model
    super_cell = model.repeat([2, 2, 2])

    #create a radial model which cuts out atoms around a representative of atom(0) in the middle of the suepr cell
    radial_model = cutout_sphere(super_cell, number_atoms_in_model, radius)
    positions_2 = radial_model.get_positions()
    number_atoms_in_radial_model = len(positions_2)
    # Create an array with the all the distances from atom the variable 'position'
    distances = radial_model.get_distances(position, range(number_atoms_in_radial_model), mic=True, vector=False)

    # plot these values as a histogram

    y, binEdges = np.histogram(distances, bins=number_atoms_in_radial_model)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-')
    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('g(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Radial Distribution Function', fontsize=15)
    pl.show()
    plt.show()
