def radial_distribution_function(model, radius, position):
    '''
    Returns a plot of the distribution of the distance between atoms up to given radius
    Plot is currently a frequency vs distance. Current usage is for periodic solids

    Parameters:

    model: Atoms objects
        The model from which the RDF is to be plotted
    radius: float
        The distance around the atom of interest to plot the RDF
    position: integer
        The atom index of interest (i.e. centre of RDF)

    Returns:

    distances: List of floats
        An sorted list of all lengths of bonds between the central atom and others in the desired radius

    TODO:   - Decide what this offers over distance_distribution_function, and make this clear.
                - It might be that we can just integrate the supercell creating into ddf.
            - Gaussian over the histogram? This is a plotting aspect
    '''

    from carmm.analyse.bonds import get_sorted_distances

    # create a super cell of the model to ensure we get all interactions within radius
    # TODO: This needs to be more rigorous e.g. check if cell is big enough
    super_cell = model.repeat([2, 2, 2])

    distances_all = get_sorted_distances(super_cell, position)
    distances = [ distance for distance in distances_all if distance < radius ]

    return sorted(distances)

def average_distribution_function(trajectory, samples=10, bin_sampling=0.1, plot=False):
    '''
    Plots the average distribution function of the last N steps of an MD trajectory
    TODO: Is this radial or distance? Should it be an option which to use?
    
    Parameters:
    
    trajectory: List of Atoms objects
        The pathway from which the ensemble RDF is to be plotted
    samples: Integer
        The number of samples to include in the ensemble, starting from the final image of the trajectory.

    Returns:

    Mean: List of floats
        The

    TODO: -OB: looks pretty complicated with all the stored variables, can this be done in a loop?
               AL: done? We need a regression test though - and comparison against the standalone script.
               All the different functions have a plot at the end, create a function that plots and can just be added

    '''

    from matplotlib import pyplot as plt
    plt.cla()
    from math import ceil
    import numpy as np
    #######
    from carmm.analyse.bonds import get_sorted_distances
    
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
        # Calculate distribution function
        snapshots_d.append(get_sorted_distances(snapshots[i]))
        # Sort data. Is this needed?
        snapshots_ds.append(sorted(snapshots_d[i]))
        # Create plot data
        y, binEdges = np.histogram(snapshots_ds[i], bins=ceil(max(snapshots_ds[i]) / bin_sampling), density=True)
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        #colors=['red','blue'] # Used for testing
        #pl.plot(bincenters, y, '-', color=colors[i], label=str(i))
        plt.plot(bincenters, y, '-', color='#808080')

    # Plot mean. This need double checking with something known e.g. H2O
    mean = [item for atoms_object in snapshots_ds for item in atoms_object]
    y, binEdges = np.histogram(mean, bins=ceil(max(mean) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    plt.plot(bincenters, y, '-', label='Mean', color="#000000")

    if plot:
        plt.xlabel('r/Å', fontsize=15)
        plt.ylabel('g(r)', fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title('Last '+str(samples)+' Snapshots', fontsize=15)
        plt.legend()
        plt.show()

    # Something should be returned.
    return mean

def plot_distribution_function(data, bins=None, bin_sampling=0.1, title=None):
    '''
    Generic plotter

    Parameters:

    data: List of floats
        The values to be plotted on the histogram
    bins: Int
        Number of separate bins on the histogram
    bin_sampling: Float
        The separation between bin boundaries
    title: String
        Name to place at the top of the plot

    TODO: Graphs need to start at zero by default
          We should edit the aces object, not the plot itself (see graphs.py, specifically mulliken related)
    '''
    from matplotlib import pyplot as plt
    from numpy import histogram

    # Calculate bins if not already defined
    if not bins:
        from math import ceil
        bins = ceil(max(data) / bin_sampling)

    y, binEdges = histogram(data, bins=bins)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    plt.plot(bincenters, y, '-')
    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('g(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(title, fontsize=15)

    # This should return the axis object, not enable it.
    #plt.show()
    return plt

