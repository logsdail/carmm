def get_graph_colour(choice=0):
    colours = ['red', 'blue', 'green', 'yellow', 'orange', 'indigo', 'violet']
    return colours[choice]

def get_graph_linetype(choice=0):
    line_types = ['solid', 'dashed', 'dashdot', 'dotted']
    return line_types[choice]

def set_graph_axes_mulliken(plt, x, y, homo, xlabel='$\epsilon$ (eV)', ylabel='Density of States (1/eV)'):
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