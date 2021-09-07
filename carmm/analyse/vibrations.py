import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.io import read


def vib_analysis(model):
    ''' Returns a graph showing displacemet of bonds/atoms in a trajectory.

    Parameters:
        model: Atoms object
               e.g trajectory file to calculate bond displacement

    TODO: Plotting should be a different function i.e plot_vibrations()
          Example needed in examples folder
    '''
    atot = Atoms.get_chemical_symbols(self=read(model))
    # print(len(atot))

    traj = Trajectory(model)

    for i in range(len(atot) - 1):
        for j in range(i + 1, len(atot)):
            distances = []
            for atoms in traj:
                dist = atoms.get_distances(i, j)
                distances.append(float(dist))
            dist_list = distances
            print(dist_list)
            x = range(len(dist_list))
            plt.plot(x, dist_list)
            plt.xlabel("time step")
            plt.ylabel('displacement')
    plt.show()