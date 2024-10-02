def get_graph_colour(choice=0):
    '''
    Return some nice graph colours

    Parameters:

    choice: Integer
        The colour choice to return from the colours array
    '''
    colours = ['red', 'blue', 'green', 'yellow', 'orange', 'indigo', 'violet']
    return colours[choice]

def get_graph_linetype(choice=0):
    '''
    Return linetypes that match the matplotlib options

    Parameters:

    choice: Integer
        The linetype choice to return from the array of options

    '''
    line_types = ['solid', 'dashed', 'dashdot', 'dotted']
    return line_types[choice]

def set_graph_axes_mulliken(axis, x, y, homo, xlabel='$\\epsilon$ (eV)', ylabel='Density of States (1/eV)'):
    '''
    Function to automate setting up the plot axes for a DOS

    Parameters:

    plt: Matplotlib plot object or Axes
        This has data predrawn from the x- and y-axis data also provided
    x: list of Floats
        The x-axis data object (i.e. range of x-axis values)
    y: list of Floats
        The y-axis data objects, corresponding to the plotted data
    homo: Float
        Reference level for the HOMO, which is the middle point of the x-axis
    xlabel: String
        Label for x-axis
    ylabel: String
        Label for y-axis
    '''

    # Determine if this is a plot instance, and if so extract axes
    from matplotlib.axes import Axes
    if(not isinstance(axis,Axes)):
        axis = axis.gca()
    # Set all the axes limits and labels
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
        axis.axhline(y=0, xmin=min(x), xmax=max(x), color='black', lw=2)
    axis.set_ylim(ymin, ymax)
    axis.set_xlim(xmin, xmax)
    axis.set_yticks([])
    # Plot HOMO line
    axis.axvline(x=homo, ymin=ymin, ymax=ymax, color='black', lw=2, ls='--')
    # Label axes
    axis.set_ylabel(ylabel)
    axis.set_xlabel(xlabel)

def load_xyz_data_from_csv(fname):
    '''
    Read in x, y and z-axis data from a CSV file for plotting (3D)

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
    Construct a suitable axes for displating a 3D heatmap

    Parameters:

    plt: Matplotlib structure
        Incoming plot containing data and for which the axes need adjustment
    x: List of Floats
        X-axis data
    y: List of Floats
        Y-axis data
    '''
    # Aesthetic values
    plt.colorbar()
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))


def plot_energy_profile(data, x_labels, **kwargs):
    '''
    This function is used for plotting solid+dashed line style reaction energy profiles, which often require manual
    tweaking using other software.

    Parameters:
        data: dict
            A dictionary of key-value pairs, where the keys are used as the title of the data series, and the values are
            lists of floats to be plotted on the y-axis.
        x_labels: list
            A list of strings containing x-axis labels - intermediates, reaction steps etc.
        **kwargs: dict
            %% These variables should be optional in the function call. Then the defaults are clear.
            This is passthrough for all other variables. The following are pulled from
            if they exist (otherwise defaults are used):
            font_size; figsize; colours; linestyles; linestyle; x_labels_rotation; x_axis_title; y_axis_title; legend_xy_offset


    Returns:
        matplotlib.pyplot
    '''

    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors
    import numpy as np

    # Endure that data is uniform for all datasets
    # TODO: allow x-y value pairs in situations where e.g. some intermediates are skipped
    for series in data:
        energies = data[series]
        assert len(energies) == len(x_labels), f"Length of {series} values: {energies} " \
                                               f"does not match the length of x_labels: " \
                                               f"{len(energies)} =/= {len(x_labels)} \n" \
                        f"This function currently only supports datasets representing an identical reaction pathway."


    # Define keyword arguments
    font_size = kwargs.get("font_size", 14)
    figsize = kwargs.get("figsize", (10, 6))
    colours = kwargs.get("colours", list(mcolors.TABLEAU_COLORS.values())[:len(data)])
    linestyles = kwargs.get("linestyles", [])
    linestyle = kwargs.get("linestyle", "-")
    x_labels_rotation = kwargs.get("x_labels_rotation", 90)
    x_axis_title = kwargs.get("x_axis_title", "{x}")
    y_axis_title = kwargs.get("y_axis_title", "Δ$E$ /eV")
    legend_xy_offset = kwargs.get("legend_xy_offset", (1, 1))

    # Create a linear plot
    fig, ax = plt.subplots(figsize=figsize)  # Set the figure size

    # Plot horizontal lines for each point and dashed lines connecting them
    for x, y_tuple in enumerate(zip(*[data[surface] for surface in data])):
        labels = [key for key in data]

        for _, y in enumerate(y_tuple):
            if linestyles:
                linestyle = linestyles[_]

            ax.hlines(y, x - 0.2, x + 0.2,
                      colors=colours[_],
                      linestyle=linestyle,
                      linewidth=2,
                      label=labels[_] if x == 0 else '')

        # Connect the points with dashed lines
        if x < len(x_labels) - 1:
            for _, y in enumerate(y_tuple):
                ax.plot([x + 0.2, x + 1 - 0.2], [y, data[labels[_]][x + 1]],
                        linestyle='--',
                        color=colours[_],
                        linewidth=1,
                        alpha=0.7)

    # Replace the existing lines for setting x-axis ticks and labels with the following:
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=x_labels_rotation, ha="center", fontsize=font_size)

    # Set y-axis font size
    ax.yaxis.set_tick_params(labelsize=font_size)

    # Set chart title and labels with increased font size
    ax.set_xlabel(x_axis_title, fontsize=font_size)
    ax.set_ylabel(y_axis_title, fontsize=font_size)

    # Add a legend with increased font size
    ax.legend(fontsize=font_size, loc='upper left', bbox_to_anchor=legend_xy_offset)

    # Add a grid
    ax.grid(axis='y', linestyle='-', alpha=0.5)

    # Remove the top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Show the chart
    plt.tight_layout()

    return plt


