'''Example of script to analyse vibrational trajectories
from a vibrational calculation using ASE'''

def test_vib_analysis():
    from carmm.analyse.vibrations import vib_analysis, plot_vibrations

    # Get H2O vib data calculated using ASE (EMT)
    file = 'data/H2O_vib/vib.1.traj'

    vib_data = vib_analysis(file)
    plot_vibrations(vib_data)

test_vib_analysis()