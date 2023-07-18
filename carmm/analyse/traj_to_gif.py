from ase.io import read
from ase.visualize import view
import os
from carmm.analyse.povray_render import povray_render, atom_sub


def traj_to_gif(filename, automatic=False, generic_projection_settings=None, povray_settings=None, frames_per_second=30,
                pause_time=0.5, atom_subs=None, convert_flags=None, keep_temp_files=True, **kwargs):
    """
    A function which takes a .traj file, visualises it in povray with your desired settings and outputs a .gif file.
    The ImageMagick linux suite is required for this function.
    :param filename: (Str) Full name/directory of the .traj file to convert
    :param automatic: (Bool) If True, automatically renders images using given settings. If False, opens ASE GUI to
    to allow the user to manually render the images (**FOLLOW THE GIVEN INSTRUCTIONS**)
    :param generic_projection_settings: (Dict) Settings used by PlottingVariables
    (see https://gitlab.com/ase/ase/-/blob/master/ase/io/utils.py PlottingVariables/__init__ for settings options)
    :param povray_settings: (Dict) Settings used by Povray for visualisation
    (see https://gitlab.com/ase/ase/-/blob/master/ase/io/pov.py POVRAY/__init__ for settings options)
    :param frames_per_second: (Float) Speed of switching images (Default is 30 fps)
    :param pause_time: (Float) Time to pause on the first and last images (Default is 0.5 seconds)
    :param atom_subs: (List of lists of strs) Pairs of atomic symbols with the first being changed to the second in
    all images for clearer visualisation
    :param convert_flags: (Dict of strs) Flags and corresponding parameters for the ImageMagick 'convert' function.
    For default parameters, leave out this parameter. Default options are '-verbose -dispose previous -loop 0'
    :param keep_temp_files: (Bool) If False, will delete the created .png files after use. If True, will keep the .png files
    and produce a .traj file if any atoms were substituted
    :return: A .gif file of the .traj file, visualised in povray
    """

    # Retrieve the file name and file extension from (potentially) a full directory path
    if filename is None:
        filename = 'atoms.traj'
    file = filename.split('/')[-1]
    file, ext = file.split('.')

    if ext != 'traj':
        raise RuntimeError('Function only supports .traj files')

    atoms = read(filename + '@:')
    steps = len(atoms)

    # Generate the list of povray image filenames
    digits = len(str(steps - 1))
    indices = [f'%0{digits}d' % i for i in range(steps)]
    filenames = [f'{file}.{index}.png' for index in indices]

    if automatic:
        for frame in range(steps):
            frame_atoms = atoms[frame]
            povray_render(frame_atoms, output=f'{file}.{indices[frame]}', view=False, atom_subs=atom_subs,
                      generic_projection_settings=generic_projection_settings,
                      povray_settings=povray_settings)
    else:
        if atom_subs is not None:
            for frame in range(steps):
                frame_atoms = atoms[frame]
                atoms[frame] = atom_subs(frame_atoms, atom_subs)
            if keep_temp_files:
                atoms.write(f'{file}_povray.traj')
        print(f'***Crucial Steps***\n'
              f'1. In ASE GUI, navigate to Tools -> Render Scene\n'
              f'2. Change "Output basename" to {file}\n'
              f'3. Select "Render all frames"\n'
              f'4. Deselect "Show output window"\n'
              f'5. Change any other settings (e.g. Atomic texture set) as desired')
        view(atoms)
        input('***Press Enter to continue once Povray is finished visualising...***\n')

    gifmaker(steps, file, filenames, frames_per_second, pause_time, convert_flags,
             keep_temp_files)

    print("Happy cooking!")

    # For testing purposes
    return file, ext, steps, atoms, filenames


def gifmaker(steps, file, filenames, frames_per_second, pause_time, convert_flags,
             keep_temp_files):

    # Get the frames_per_second into a format that ImageMagick accepts for the delay flag
    if frames_per_second <= 1:
        delay = str(1 / frames_per_second)
    else:
        delay = f'1x{frames_per_second}'

    # Add in extra first/last frames at the start/end to give a pause of the desired length
    if pause_time is not None:
        pause_frames = int(pause_time * frames_per_second)
        count = 1  # Already one frame in the list
        while count in range(pause_frames):
            filenames.insert(0, filenames[0])
            filenames.insert(-1, filenames[-1])
            count += 1

    # Convert flags
    # Default values
    if convert_flags is None:
        convert_flags = {'-verbose': '', '-dispose': 'previous', '-loop': '0'}
    # Concatenate string for command
    convert_options = ''
    for key, value in convert_flags.items():
        convert_options += f'{key} {value} '

    # Execute the ImageMagick convert command in the terminal
    try:
        command = (f'convert {convert_options} -delay {delay} %s {file}.gif' % ' '.join(filenames))
        os.system(command)
    except:
        print('FileNotFoundError: Files do not exist')

    # Delete the povray image files if requested
    if not keep_temp_files:
        for index in filenames:
            os.system(f'rm {filenames[index]}')

    # For testing purposes
    return filenames, delay, convert_options
