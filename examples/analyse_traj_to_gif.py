'''
Example script showing how to convert an example trajectory file (.traj)
into a gif visualised in povray.
'''

def test_traj_to_gif():
    from carmm.analyse.traj_to_gif import traj_to_gif, gifmaker, atom_sub
    from ase.io import read

    # Call the overall function like this, the rest of the script just tests the functions individually
    # traj_to_gif('data/NH3-H3O_traj/nh3-h3o.traj@:', frames_per_second=10,
    #             pause_time=1, atom_subs=[['N', 'C']], keep_temp_files=True)

    file = 'data/NH3-H3O_traj/nh3-h3o.traj@:'

    atoms = read(file)
    frame_atoms_list = atom_sub(atoms, [['N', 'C']], 41, 'nh3-h3o', 'traj')
    assert len(frame_atoms_list) == 41
    assert frame_atoms_list[11].symbols[4] == 'C'  # Random frame

    filenames, delay = gifmaker(41, 'nh3-h3o', 'traj', frames_per_second=10, pause_time=1, keep_png_files=False)
    assert delay == '1x10'
    assert len(filenames) == 59
