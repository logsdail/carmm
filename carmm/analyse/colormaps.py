def color_bar(width, height, bottom, vmin, vmax, orientation, label, cmapcolour):
    """Basic Colour Bar
    figsize= y axis how tall
    bottom = <0.8 - size of bar (cannot be bigger than top)
    vmin/vmax = x axis values
    orientation = horizontal,vertical
    label = units
                                    """

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import cm

    fig, ax = plt.subplots(figsize=(width, height))
    fig.subplots_adjust(bottom=bottom)
    colour = cm.get_cmap(cmapcolour, 1000)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colour),
                cax=ax, orientation=orientation, label=label)

    return fig

