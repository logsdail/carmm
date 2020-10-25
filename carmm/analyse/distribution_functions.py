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

def average_distribution_function(trajectory):

   '''plots the average distribution function of the last 10 steps of an MD trajectory

   TODO: -OB: looks pretty complicated with all the stored variables, can this be done in a loop?
               All the different functions have a plot at the end, create a function that plots and can just be added '''



    from ase.io import read
    from matplotlib import pyplot as plt
    from math import ceil
    import numpy as np
    import pylab as pl
    import distance_distribution_function as ddf

    # Variable Naming:
    # [a-l] - Snapshot
    # [a-l]d - Distance distribution of the snapshot
    # [a-l]ds - Distance distribution array sorted from lowest-highest

    bin_sampling = 0.1

    a = read(trajectory, -1)
    b = read(trajectory, -2)
    c = read(trajectory, -3)
    d = read(trajectory, -4)
    e = read(trajectory, -5)
    f = read(trajectory, -6)
    g = read(trajectory, -7)
    h = read(trajectory, -8)
    l = read(trajectory, -9)
    j = read(trajectory, -10)

    del a.constraints
    del b.constraints
    del c.constraints
    del d.constraints
    del e.constraints
    del f.constraints
    del g.constraints
    del h.constraints
    del l.constraints
    del j.constraints

    ad = ddf(a)
    bd = ddf(b)
    cd = ddf(c)
    dd = ddf(d)
    ed = ddf(e)
    fd = ddf(f)
    gd = ddf(g)
    hd = ddf(h)
    ld = ddf(l)
    jd = ddf(j)

    ads = sorted(ad)
    bds = sorted(bd)
    cds = sorted(cd)
    dds = sorted(dd)
    eds = sorted(ed)
    fds = sorted(fd)
    gds = sorted(gd)
    hds = sorted(hd)
    lds = sorted(ld)
    jds = sorted(jd)

    mean = ads + bds + cds + dds + eds + fds + gds + hds + lds + jds
    np.true_divide(mean, 10)

    y, binEdges = np.histogram(ads, bins=ceil(max(ads) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(bds, bins=ceil(max(bds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(cds, bins=ceil(max(cds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(dds, bins=ceil(max(dds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(eds, bins=ceil(max(eds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(fds, bins=ceil(max(fds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(gds, bins=ceil(max(gds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(hds, bins=ceil(max(hds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(lds, bins=ceil(max(lds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(jds, bins=ceil(max(jds) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', color="#808080")
    y, binEdges = np.histogram(mean, bins=ceil(max(mean) / bin_sampling), density=True)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    pl.plot(bincenters, y, '-', label='Mean', color="#000000")

    plt.xlabel('r/Å', fontsize=15)
    plt.ylabel('g(r)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Last 10 Snapshots', fontsize=15)
    pl.legend()
    pl.show()
    plt.show()

