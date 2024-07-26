"""
Test script for traj_to_gif and gifmaker functions.
Usage example in comment block, showing how to convert an example trajectory file (.traj)
into a gif visualised in povray.
"""

import os

def test_traj_to_gif():
    from carmm.utils.traj_to_gif import traj_to_gif, gifmaker

    # Usage example:
    # Call the overall function like this, the rest of the script just tests the functions individually
    # traj_to_gif('data/NH3-H3O_traj/nh3-h3o.traj', frames_per_second=10,
    #             pause_time=1, atom_subs=[['N', 'C'], ['O', 'N']], keep_temp_files=True)

    file = 'data/NH3-H3O_traj/nh3-h3o.traj'
    file, ext, steps, atoms, filenames = traj_to_gif(file, automatic=True, pause_time=None,
                                                     atom_subs=[['N', 'C'], ['O', 'N']])

    assert file == 'nh3-h3o'
    assert ext == 'traj'
    assert steps == 41
    assert len(atoms[17]) == 8  # Random frame
    assert atoms[11].symbols[4] == 'C'  # Random frame
    assert atoms[38].symbols[0] == 'N'  # Random frame
    assert filenames[7] == 'nh3-h3o.07.png'  # Random frame

    indices = [f'%02d' % i for i in range(steps)]
    filenames_gif, duration, gif_options = gifmaker('nh3-h3o', filenames=filenames, frames_per_second=10,
                                                    pause_time=1, gif_options=None, indices=indices,
                                                    keep_temp_files=True)
    assert len(filenames_gif) == 41
    assert duration == [1000] + [100]*(steps-2) + [1000]
    assert gif_options == {
        'save_all': True,
        'optimize': False,
        'loop': 0,
    }

    gif_options = {
        'optimize': True,
        'loop': 2,
        'comment': 'A random other option'
    }
    filenames = [f'{file}.{index}.png' for index in indices]  # Redefine 'filenames' for new pause_time
    filenames_gif, duration, gif_options = gifmaker('nh3-h3o', filenames=filenames, frames_per_second=0.1,
                                                     pause_time=20, gif_options=gif_options, indices=indices,
                                                     keep_temp_files=False)
    assert len(filenames_gif) == 41
    assert duration == [20000] + [10000]*(steps-2) + [20000]
    assert gif_options == {
        'optimize': True,
        'loop': 2,
        'comment': 'A random other option',
        'save_all': True,
    }


test_traj_to_gif()

os.system('rm nh3-h3o.gif')
