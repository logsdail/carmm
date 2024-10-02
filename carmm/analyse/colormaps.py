def color_bar(width, height, bottom, vmin, vmax, orientation, label, cmapcolour, resolution):
    """Basic Colour Bar\
    figsize= y axis how tall\  
    bottom = <0.8 - size of bar (cannot be bigger than top)\  
    vmin/vmax = x axis values\  
    orientation = horizontal,vertical\  
    label = units
    
    cmapcolour = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma', 'cividis']),
         ('Sequential', [
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
         ('Sequential (2)', [
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper']),
         ('Diverging', [
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
         ('Cyclic', ['twilight', 'twilight_shifted', 'hsv']),
         ('Qualitative', [
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c']),
         ('Miscellaneous', [
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
            'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral',
            'gist_ncar'])]
    resolution = the level of segmentation on the map
                                    """

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import cm

    fig, ax = plt.subplots(figsize=(width, height))
    fig.subplots_adjust(bottom=bottom)
    colour = cm.get_cmap(cmapcolour, resolution)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colour),
                 cax=ax, orientation=orientation, label=label)

    return fig
