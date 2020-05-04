#!/usr/bin/env python3

'''
This is a more complicated example and QA test for extracting and plotting Mulliken data from FHI-aims
Note that this routine plots a graph with both spin-up and spin-down data, but also extracts
information specific to atoms and angular momenta of interest, and plots all the data in separate graphs

This is useful when trying to understand the electronic structure of your model
'''

def test_analyse_mulliken_atomic_comparison():

    import matplotlib.pyplot as plt
    from carmm.analyse.mulliken import parse_mulliken_file
    from carmm.analyse.graphs import get_graph_linetype, get_graph_colour, set_graph_axes_mulliken
    from carmm.analyse.bonds import get_indices_of_elements

    output_files = ['data/CO/co_light.log', 'data/Fe/fe_light.log', 'data/Fe-CO/fe-co_light.log']
    mulliken_files = ['data/CO/Mulliken.out', 'data/Fe/Mulliken.out', 'data/Fe-CO/Mulliken.out']

    # Store the axes for different graphs - this is a matplotlib specific aspect.
    # Instead of drawing in data on to the plt object, here we plot on to the axes objects.
    fig, axes = plt.subplots(len(output_files),1)

    for i in range(len(output_files)):
        # Read in atoms information
        output_file = output_files[i]
        from ase.io import read
        atoms = read(output_file)

        # Read in Mulliken data from file
        mulliken_file = mulliken_files[i]
        mulliken_data = parse_mulliken_file(mulliken_file)

        # Collect all the density of states data to plot
        x, all_data = mulliken_data.get_all_plot_data()

        # Collect the indices for each element we are interested in
        fe_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'fe')
        c_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'c')
        o_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'o')

        # Collect the density of states data to plot as a function of atomic label
        x, fe = mulliken_data.get_plot_data(fe_indices, range(mulliken_data.get_nspin()),
                                            range(mulliken_data.get_nkpts()), 'd')
        x, c = mulliken_data.get_plot_data(c_indices, range(mulliken_data.get_nspin()),
                                            range(mulliken_data.get_nkpts()), 'sp')
        x, o = mulliken_data.get_plot_data(o_indices, range(mulliken_data.get_nspin()),
                                            range(mulliken_data.get_nkpts()), 'sp')

        # The returned array has either one entries if spin-collinear, or two entries for spin up and down.
        for spin in range(mulliken_data.get_nspin()):
            if spin == 0:
                # Here we are filling between lines, rather than plotting a single line.
                # Hence we provide two y-values to plot between.
                axes[i].fill_between(x, [ 0 * len(x) ], fe[spin], lw=0, facecolor=get_graph_colour(0), label='Fe', interpolate=True)
                axes[i].fill_between(x, fe[spin], fe[spin]+c[spin], lw=0, facecolor=get_graph_colour(1), label='C', interpolate=True)
                axes[i].fill_between(x, fe[spin]+c[spin], fe[spin]+c[spin]+o[spin], lw=0, facecolor=get_graph_colour(2), label='O', interpolate=True)
                axes[i].plot(x, all_data[spin], lw=1, color='black', ls=get_graph_linetype())
            else: # (spin == 1)
                # Note the slightly different definition of the y-axis data here
                # There is a chance that the element doesn't exist in the system, which returns a DOS of zero,
                # and you cannot take the negative of a list of zeros (it transpires)
                # Instead, the approach given just makes each value in the list its negative value separately
                axes[i].fill_between(x, [0 * len(x)], [-i for i in fe[spin]], lw=0, facecolor=get_graph_colour(0), interpolate=True)
                axes[i].fill_between(x, [-i for i in fe[spin]], [-i for i in fe[spin]+c[spin]], lw=0, facecolor=get_graph_colour(1), interpolate=True)
                axes[i].fill_between(x, [-i for i in fe[spin]+c[spin]], [-i for i in fe[spin]+c[spin]+o[spin]], lw=0, facecolor=get_graph_colour(2), interpolate=True)
                axes[i].plot(x, -(all_data[spin]), lw=1, color='black', ls=get_graph_linetype())

        # Work to rescale axes. Extracts the maximum y-value
        set_graph_axes_mulliken(axes[i], x, all_data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

        # Add a legend
        axes[i].legend()

    #### Shoe-horn an assertion test in on final HOMO ####
    assert(mulliken_data.get_homo() == -4.22285)
    ########

    # Display the graphs
    # plt.show()

# Run the example
test_analyse_mulliken_atomic_comparison()