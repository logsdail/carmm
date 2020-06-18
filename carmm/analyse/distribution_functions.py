def distance_distribution_function(model):
    '''Returns a plot of the distribution of the distance between all atoms
    plot is currently a frequency vs distance. Current usage is for periodic solids
   TODO: - currently plot is a histogram need to change to be a gaussian - separate the bins?
           AJL: Could the plot be at the taste of the user, and the function just returns the data?
           You could perhaps have a default x-axis or spacing of bins, and then this could be overwritten by user?
         - Only calculates RDF with respect the first atom - needs to be generalised:
            - what about for any other atom of interest?
            - what about averaging over all atoms, as per standard EXAFS?
              OB: amended distance_distribution to average over all atoms

        '''
    from matplotlib import pyplot as plt
    # get all distances in the model
    distances = model.get_all_distances(mic=True, vector=False)
    # current distance variable is an array of coordinates for each atom, average position for each atom:
    a = (sum(distances) / 3)
    # get the relative positions of atoms to eachother
    norm = a - min(a)

    # plot these values as a histogram
    plt.figure()
    plt.hist(norm, bins=len(distances), density=False, histtype='step', align='mid', stacked=True)
    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('n(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Distribution Function', fontsize=15)
    plt.show()


def radial_distribution_function(model, radius, position):
    '''Returns a plot of the distribution of the distance between each atom from atom_0.
    plot is currently a frequency vs distance. Current usage is for periodic solids
    This script will create a radius around a given model and calculate the distance of this new model
    TODO:  - Gaussian over the histogram
            AJL: What's the difference between this and difference_distribution_function?
    '''
    from ase.io import read
    from matplotlib import pyplot as plt
    from carmm.build.cutout import cutout_sphere

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
    plt.figure()
    plt.hist(distances ,bins=number_atoms_in_radial_model, density=False, color='b', histtype='step', align='mid', stacked=True)
    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('g(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Radial Distribution Function', fontsize=15)

    #plot data as gaussian function

    plt.show()