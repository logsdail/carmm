def get_graph_colour(choice=0):
    '''
    Description

    Parameters:

    choice: Integer
        The colour choice to return from the colours array
    '''
    colours = ['red', 'blue', 'green', 'yellow', 'orange', 'indigo', 'violet']
    return colours[choice]

def get_graph_linetype(choice=0):
    '''
    Description

    Parameters:

    choice: Integer
        The linetype choice to return from the array of options

    '''
    line_types = ['solid', 'dashed', 'dashdot', 'dotted']
    return line_types[choice]

def set_graph_axes_mulliken(plt, x, y, homo, xlabel='$\epsilon$ (eV)', ylabel='Density of States (1/eV)'):
    '''
    Description

    Parameters:

    plt:

    x:

    y:

    homo:

    xlabel:

    ylabel:

    '''
    ymax = max(map(max, y))*1.1
    ymin = 0
    if homo == 0.0:
        xmax = 10
        xmin = -10
    else:
        xmax = 10 #max(x)
        xmin = -15 #min(x)
    if len(y) > 1:
        ymin = -ymax
        # Add zero line for spin-polarised systems
        plt.axhline(y=0, xmin=min(x), xmax=max(x), color='black', lw=2)
    plt.ylim(ymin, ymax)
    plt.xlim(xmin, xmax)
    plt.yticks([])
    # Plot HOMO line
    plt.axvline(x=homo, ymin=ymin, ymax=ymax, color='black', lw=2, ls='--')
    # Label axes
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

def load_xyz_data_from_csv(fname):
    '''
    Description

    Parameters:

    fname: String
        filename for the CSV that is to be read in
    '''
    import numpy as np
    # Load data from csv
    dat = np.genfromtxt(fname, delimiter=', ',skip_header=0)
    X_dat = dat[:,0]
    Y_dat = dat[:,1]
    Z_dat = dat[:,2]

    # Convert from pandas dataframes to numpy arrays
    X, Y, Z, = np.array([]), np.array([]), np.array([])
    for i in range(len(X_dat)):
        X = np.append(X,X_dat[i])
        Y = np.append(Y,Y_dat[i])
        Z = np.append(Z,Z_dat[i])

    return X, Y, Z

def set_graph_axes_heatmap(plt, x, y):
    '''
    Description

    Parameters:

    plt:

    x:

    y:

    '''
    # Aesthetic values
    plt.colorbar()
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))