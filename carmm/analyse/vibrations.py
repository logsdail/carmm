import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.io import read


def vib_analysis(model):
    ''' Returns a list displacemet of bonds/atoms in a trajectory.
    Parameters:
        model: Atoms object
               e.g trajectory file to calculate bond displacement
    TODO: - Resolve for periodic systems - functionally currently doesn't work


    '''
    atot = Atoms.get_chemical_symbols(self=read(model))

    traj = Trajectory(model)

    for i in range(len(atot) - 1):
        for j in range(i + 1, len(atot)):
            distances = []
            for atoms in traj:
                dist = atoms.get_distances(i, j)
                distances.append(float(dist))
            dist_list = distances
            return dist_list


class plot_vibration_data:
    ''' Returns a graph showing displacement of bonds/atoms in a vibration trajectory from ASE.
    #     Parameters:
    #         x_axis: length of data returned from vib_analysis()
              y_xis: data returned by vib_analysis()
             title: title for plot

             TODO: would be nice to able to plot atomic/elemental information on plot i.e which atoms are being displaced
    '''


    def __init__(self,x_axis, y_axis, title):
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.title = title

    def plot_vib(self):
        plt.plot(self.x_axis, self.y_axis)
        plt.title(self.title)
        plt.show()
