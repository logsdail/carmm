def distance_distribution_function(model, bin_sampling, plot=False):
    '''Returns a plot of the distribution of the distance between all atoms
    plot is currently a frequency vs distance. Current usage is for periodic solids
   TODO:
         - This needs a list of what the inputs/outputs are. The documentation isn't good enough!
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
    if plot:
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

    return individual_lengths

def radial_distribution_function(model, radius, position, plot=False):
    '''Returns a plot of the distribution of the distance between each atom from atom_0.
    plot is currently a frequency vs distance. Current usage is for periodic solids
    This script will create a radius around a given model and calculate the distance of this new model
    TODO:   - This needs a list of inputs/outputs, as otherwise it is really tough to work with.
            - Gaussian over the histogram
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
    if plot:
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

    return distances

def average_distribution_function(trajectory, samples=10, bin_sampling=0.1, plot=False):
    '''
    Plots the average distribution function of the last N steps of an MD trajectory
    TODO: Is this radial or distance? Should it be an option which to use?
    
    Parameters:
    
    trajectory: List of Atoms objects
        The pathway from which the ensemble RDF is to be plotted
    samples: Integer
        The number of samples to include in the ensemble, starting from the final image of the trajectory.

    TODO: -OB: looks pretty complicated with all the stored variables, can this be done in a loop?
               AL: done? We need a regression test though - and comparison against the standalone script.
               All the different functions have a plot at the end, create a function that plots and can just be added '''

    from ase.io import read
    from matplotlib import pyplot as plt
    from math import ceil
    import numpy as np
    import pylab as pl
    
    # Suggested update so we can use loops:
    # snapshots[0-n] - snapshots
    # snapshots_d[0-n] - Distance distribution of snapshots
    # snapshots_ds[0-n] - Distance snapshots sorted from lowest to highest
    
    snapshots = []
    snapshots_d = []
    snapshots_ds = []
    
    # Read in all Atoms objects to be sampled
    for i in range(samples):
        snapshots.append(trajectory[(-i)-1])
        # Remove any constraints, as these hamper analysis
        del snapshots[i].constraints
        # Caculate distribution function
        # @OB: distance_distribution_function doesn't return anything, so this setup doesn't work
        #      Was the code edited for ddf by Corey? If so this needs incorporating
        #      (Really DDF should return the data anyway, not plot it - the plotter should be in graphs.py)
        snapshots_d.append(distance_distribution_function(snapshots[i], bin_sampling))
        # Sort data. Is this needed?
        snapshots_ds.append(sorted(snapshots_d[i]))
        # Create plot data
        y, binEdges = np.histogram(snapshots_ds[i], bins=ceil(max(snapshots_ds[i]) / bin_sampling), density=True)
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        #colors=['red','blue'] # Used for testing
        #pl.plot(bincenters, y, '-', color=colors[i], label=str(i))
        pl.plot(bincenters, y, '-', color='#808080')

    # Plot mean. This need double checking with something known e.g. H2O
    mean = [item for atoms_object in snapshots_ds for item in atoms_object]
    y, binEdges = np.histogram(mean, bins=ceil(max(mean) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', label='Mean', color="#000000")

    if plot:
        plt.xlabel('r/Å', fontsize=15)
        plt.ylabel('g(r)', fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title('Last '+str(samples)+' Snapshots', fontsize=15)
        # What is the difference between pl and plt? This needs sorting.
        pl.legend()
        pl.show()
        plt.show()

    # Something should be returned.

