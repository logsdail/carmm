def test_plot_energy_profile():
    '''

    Returns:

    '''
    from carmm.analyse.graphs import plot_energy_profile
    import random

    # Provide the x-axis labels (x-axis is the reaction coordinate by default)
    # The labels can contain LaTeX-style formatting
    states = [fr"x$_{n}$" for n in range(6)]

    # Generate random data for plotting
    dummy_data_dictionary = {f"Catalyst {n}": [random.uniform(-1, 1.5) for i in range(len(states))]
                  for n in range(2)}

    plt = plot_energy_profile(data=dummy_data_dictionary, x_labels=states)

    # Assert number of plotted lines remains consistent
    assert len(plt.gca().get_lines()) / (len(states) - 1) == len(dummy_data_dictionary)

    # plt.show()


test_plot_energy_profile()