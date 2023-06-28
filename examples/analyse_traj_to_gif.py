"""
Example script showing how to convert an example trajectory file (.traj)
into a gif visualised in povray.
"""


def test_traj_to_gif():
    from carmm.analyse.traj_to_gif import traj_to_gif, gifmaker, povray_render
    from ase.io import read

    # Call the overall function like this, the rest of the script just tests the functions individually
    # traj_to_gif('data/NH3-H3O_traj/nh3-h3o.traj', frames_per_second=10,
    #             pause_time=1, atom_subs=[['N', 'C'], ['O', 'N']], keep_temp_files=True)

    file = 'data/NH3-H3O_traj/nh3-h3o.traj'
    file, ext, steps = traj_to_gif(file, test=True)
    assert file == 'nh3-h3o'
    assert ext == 'traj'
    assert steps == 41

    atoms = read('data/NH3-H3O_traj/nh3-h3o.traj@:')
    frame_atoms_list = povray_render(atoms, [['N', 'C'], ['O', 'N']], 41, 'nh3-h3o', 'traj', test=True)
    assert len(frame_atoms_list) == 41
    assert frame_atoms_list[11].symbols[4] == 'C'  # Random frame
    assert frame_atoms_list[38].symbols[0] == 'N'  # Random frame

    fps = 10
    filenames, delay = gifmaker(41, 'nh3-h3o', 'traj', frames_per_second=fps, pause_time=1, convert_flags=None,
                                keep_png_files=False, test=True)
    assert len(filenames) == 59
    assert delay == '1x10'

    fps = 0.1
    filenames, delay = gifmaker(41, 'nh3-h3o', 'traj', frames_per_second=fps, pause_time=1, convert_flags=None,
                                keep_png_files=False, test=True)
    assert delay == '10.0'


test_traj_to_gif()
