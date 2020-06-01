
def distance_distributon_function(model):
    '''Returns a plot of the distribution of the distance between each atom from atom_0.
    plot is currently a frequency vs distance. Current usage is for periodic solids
    TODO: - normalise distribution, atoms at the opposite end of the cell are actually close together - need to account for this.
            AJL: I think get_distances(mic=True) will help with this? Maybe from the Atoms object (https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.get_distance)
            or also from the Geometry object: https://wiki.fysik.dtu.dk/ase/ase/geometry.html#ase.geometry.get_distances
          - currently plot is a histogram need to chnage to be a gaussian
            AJL: Could the plot be at the taste of the user, and the function just returns the data?
            You could perhaps have a default x-axis or spacing of bins, and then this could be overwritten by user?
    '''
    from ase.io import read
    from matplotlib import pyplot as plt
    
    #Read file or Atoms object
    if isinstance(model, str) is True:
        model = read(model)

     #Create a variable which represents the amount of atoms in the system
    positions = model.get_positions()
    A = len(positions)
    #Create an array with the all the distances from atom 0
    distances = model.get_distances(0, range(A), mic=False, vector=False)
    #plot these values as a histogram
    plt.figure()
    plt.hist(distances ,bins=A, density=False, color='b', histtype='step', align='mid', stacked=True)
    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('n(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Distribution Function', fontsize=15)
    plt.show()


def radial_distributon_function(model):
    '''Returns a plot of the distribution of the distance between each atom from atom_0.
    plot is currently a frequency vs distance. Current usage is for periodic solids
     TODO: - Ammend script to account for a radial function
            - Normalise distribution
            - Gaussians over the histrogram
            AJL: What's the difference between this and difference_distribution_function?
    '''
    from ase.io import read
    from matplotlib import pyplot as plt
   
    #Read file or Atoms object
    if isinstance(model, str) is True:
        model = read(model)

    # Create a variable which represents the amount of atoms in the system
    positions = model.get_positions()
    A = len(positions)
    # Create an array with the all the distances from atom 0
    distances = model.get_distances(0, range(A), mic=False, vector=False)
    # plot these values as a histogram
    plt.figure()
    plt.hist(distances ,bins=A, density=False, color='b', histtype='step', align='mid', stacked=True)
    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('g(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Radial Distribution Function', fontsize=15)
    plt.show()
