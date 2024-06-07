def radial_distribution_function(model, radius, position, verbose=False):
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

    '''

    from carmm.analyse.bonds import get_sorted_distances
    from math import ceil

    # Create a super cell of the model big enough to get all interactions within radius
    lattice_vectors = model.cell.cellpar()[:3]
    super_cell_repeat = [int(ceil((radius/lattice_vectors[0])*2)),
                         int(ceil((radius/lattice_vectors[1])*2)),
                         int(ceil((radius/lattice_vectors[2])*2))]
    if verbose:
        print("A, B and C vector lengths = ", lattice_vectors)
        print("Super cell expansion necessary to accommodate radius = ", super_cell_repeat)
    super_cell = model.repeat(super_cell_repeat)

    distances_all = get_sorted_distances(super_cell, position)
    distances = [ distance for distance in distances_all if distance < radius ]

    return sorted(distances)

def element_radial_distribution_function(model, radius, element,position=None, verbose=False):
    '''
    Returns a plot of the distribution of the distance between atoms of a specific element up to given radius
    Plot is currently a frequency vs distance. Current usage is for periodic solids

    Parameters:

    model: Atoms objects
        The model from which the RDF is to be plotted
    radius: float
        The distance around the atom of interest to plot the RDF
    position: integer
        The atom index of interest (i.e. centre of RDF). Default is 0

    Returns:

    distances: List of floats
        An sorted list of all lengths of bonds between the central atom and others in the desired radius

    '''

    from carmm.analyse.bonds import get_sorted_distances
    from math import ceil
    from ase.visualize import view
    #Get chemical symbol
    symbol = model.get_chemical_symbols()

    #Get indicies of element of interest
    indices = [i for i, x in enumerate(symbol) if x == element]

    #Update atoms object for only element of interest
    model = model[indices]

   #If position is not set will automatically use index 0, otherwise will used specified position
    if position is None:
        position = 0

    # Create a super cell of the model big enough to get all interactions within radius
    lattice_vectors = model.cell.cellpar()[:3]
    super_cell_repeat = [int(ceil((radius/lattice_vectors[0])*2)),
                         int(ceil((radius/lattice_vectors[1])*2)),
                         int(ceil((radius/lattice_vectors[2])*2))]
    from ase.visualize import view

    if verbose:
        print("A, B and C vector lengths = ", lattice_vectors)
        print("Super cell expansion necessary to accommodate radius = ", super_cell_repeat)
    super_cell = model.repeat(super_cell_repeat)

    distances_all = get_sorted_distances(super_cell, position)
    distances = [distance for distance in distances_all if distance < radius ]

    return sorted(distances)

def average_distribution_function(trajectory, samples=10):
    '''
    Plots the average distribution function of the last N steps of an MD trajectory
    TODO: -This is distance based - can we adapt it to also allow radial?
          - Seems like constraints are removed below but this causes error
          when no constraints are present in the system. Need to make this a condition.
    
    Parameters:
    
    trajectory: List of Atoms objects
        The pathway from which the ensemble RDF is to be plotted
    samples: Integer
        The number of samples to include in the ensemble, starting from the final image of the trajectory.

    Returns:

    all_data: List of floats
        A list containing _all_ of the distances encountered in the sampled trajectories
    snapshots: List of list of floats
        A list containing the list of floats for distances measured in each specific snapshot analysed

    '''

    from carmm.analyse.bonds import get_sorted_distances

    # snapshots[0-n] - Distance snapshots sorted from lowest to highest
    snapshots = []
    
    # Read in all Atoms objects to be sampled
    for i in range(samples):
        model = trajectory[(-i)-1].copy()
        # Remove any constraints, as these hamper analysis.
        # Can't see why this is needed? A more verbose comment would help. TODO
        del model.constraints
        # Calculate distribution function
        snapshots.append(get_sorted_distances(model))

    # Plot mean. This need double checking with something known e.g. H2O
    all_data = [item for atoms_object in snapshots for item in atoms_object]

    # Returning the mean and the snapshot data.
    return all_data, snapshots

def plot_distribution_function(data, bins=None, bin_sampling=0.1, title=None, density=False, **kwargs):
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
    **kwargs: Additional arguments
        These are passed straight through to the plotting function

    TODO: We should edit the axes object, not the plot itself (see graphs.py, specifically mulliken related)
    '''
    from matplotlib import pyplot as plt
    from numpy import histogram

    # Calculate bins if not already defined
    if not bins:
        from math import ceil
        bins = ceil(max(data) / bin_sampling)

    y, binEdges = histogram(data, bins=bins, density=density)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    plt.plot(bincenters, y, '-', **kwargs)
    plt.xlim(xmin=0)

    # Aesthetic
    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('g(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(title, fontsize=15)

    # This should return the axis object, not enable it.
    return plt

def radius_of_gyration(model):
    ''' creates a radius of gyration for input file of a molecule.
        A radius of gyration is a measurement of the distribution of atoms
        in a molecular structure with respect to it's centre of mass: 

        radius_of_gyration^2 = SUM(mass_atom(atom_position - centre_mass_position)) / mass_of_molecule

        TODO: - further tidying up
     '''
    import numpy as np

    # creating variables from calulations using ase atoms class
    mass_array = model.get_masses()
    position_array = model.get_positions()
    center_mass = model.get_center_of_mass()

    # calculates mass of molecule
    mass_of_molecule = np.sum(mass_array)

    # calculate distance from centre of mass for all atoms and multiplies them by the mass of their respective atom
    lg = len(position_array)
    m = np.zeros(shape=(lg,3))
    for i in range(lg):
        m[i] = mass_array[i] * (position_array[i] - center_mass)**2
		
    # calculates the radius of gyration
    rog2 = np.sum(m) / mass_of_molecule
    radius_gyration = np.sqrt(rog2)

    return radius_gyration


