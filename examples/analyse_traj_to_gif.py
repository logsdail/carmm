"""
Test script for traj_to_gif and gifmaker functions.
Usage example in comment block, showing how to convert an example trajectory file (.traj)
into a gif visualised in povray.
"""


def test_traj_to_gif():
    from carmm.analyse.traj_to_gif import traj_to_gif, gifmaker
    from ase.io import read

    # Usage example:
    # Call the overall function like this, the rest of the script just tests the functions individually
    # traj_to_gif('data/NH3-H3O_traj/nh3-h3o.traj', frames_per_second=10,
    #             pause_time=1, atom_subs=[['N', 'C'], ['O', 'N']], keep_temp_files=True)

    file = 'data/NH3-H3O_traj/nh3-h3o.traj'
    file, ext, steps, atoms, filenames = traj_to_gif(file, automatic=True, atom_subs=[['N', 'C'], ['O', 'N']])

    assert file == 'nh3-h3o'
    assert ext == 'traj'
    assert steps == 41
    assert len(atoms[17]) == 8  # Random frame
    assert frame_atoms_list[11].symbols[4] == 'C'  # Random frame
    assert frame_atoms_list[38].symbols[0] == 'N'  # Random frame
    assert filenames[7] == 'nh3-h3o.07.png'

    fps = 10
    filenames_gif, delay, convert_options = gifmaker(41, 'nh3-h3o', filenames=filenames, frames_per_second=fps,
                                                     pause_time=1, convert_flags=None, keep_temp_files=False)
    assert len(filenames_gif) == 59
    assert delay == '1x10'
    assert convert_options == '-verbose -dispose previous -loop 0'

    fps = 0.1
    convert_flags = {
        '-flip': '',
        '-loop': 2,
        '-autogamma': '',
    }
    filenames_gif, delay, convert_options = gifmaker(41, 'nh3-h3o', filenames=filenames, frames_per_second=fps,
                                                     pause_time=20, convert_flags=convert_flags, keep_temp_files=True)
    assert len(filenames_gif) == 43
    assert delay == '10.0'
    assert convert_options == '-flip -loop 2 -autogamma'


test_traj_to_gif()
