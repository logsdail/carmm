'''Example of script to analyse vibrational trajectories
from a vibrational calculation using ASE'''

def test_vib_analysis():
    from carmm.analyse.vibrations import vib_analysis, plot_vibration_data

    # Get H2O vib data calculated using ASE (EMT)
    file = 'data/H2O_vib/vib.1.traj'

    vib_data = vib_analysis(file)

    ## plot data using plot class

    x = range(len(vib_data))
    y = vib_data

    plot = plot_vibration_data(x,y,'example vib_analysis')
    plot.plot_vib()

test_vib_analysis()
