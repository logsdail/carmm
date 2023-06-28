from ase.io import read
from ase.visualize import view
import os


def traj_to_gif(filename, frames_per_second=30, pause_time=0.5, atom_subs=None, convert_flags=None,
                keep_temp_files=True, test=False):

    """
    A function which takes a .traj file, visualises it in povray with your desired settings and outputs a .gif file.
    When the function gives you the ase viewer, rotate to your desired view, go to Tools -> Render Scene, select the
    "Render all frames" option and deselect "Show output window".
    MAKE SURE THAT THE OUTPUT BASENAME IS THE SAME AS THE INITIAL FILENAME (e.g. for 'atoms.traj' the output basename is
    'atoms').
    The rest of the settings can be changed at your discretion.
    Then press "Render" and wait for it to finish before closing the window and pressing enter to continue.
    The ImageMagick linux suite is required for this function.
    :param filename: (str) Full name/directory of the .traj file to convert
    :param frames_per_second: (float) Speed of switching images (Default is 30 fps)
    :param pause_time: (float) Time to pause on the first and last images (Default is 0.5 seconds)
    :param atom_subs: (list of lists of strings) Pairs of atomic symbols with the first being changed to the second in
    all images for clearer visualisation
    :param convert_flags: (dict of strings) Flags and corresponding parameters for the ImageMagick 'convert' function.
    For default parameters, leave out this parameter. Default options are '-verbose -dispose previous -loop 0'
    :param keep_temp_files: (boolean) If set to False, will delete the created povray .traj and .png files after use
    :param test: (boolean) Do not change! A quick fix to stop the unittest from getting stuck waiting for a user input
    :return: A .gif file of the .traj file, visualised in povray
    """

    # Retrieve the file name and file extension from (potentially) a full directory path
    if filename is None:
        filename = 'atoms.traj'
    file = filename.split('/')[-1]
    file, ext = file.split('.')

    atoms = read(filename+'@:')
    steps = len(atoms)

    if not test:
        povray_render(atoms, steps, file, ext, atom_subs)

    gifmaker(steps, file, ext, frames_per_second, pause_time, convert_flags, keep_temp_files, test)

    print("Happy cooking!")

    # For testing purposes
    return file, ext, steps


def atom_sub(atoms, atom_subs, steps, file, ext):

    frame_atoms_list = []

    # Replace atoms of one element with another for clearer visualisation
    for frame in range(steps):
        frame_atoms = atoms[frame]
        if atom_subs is not None:
            for sub in atom_subs:
                for i in range(len(frame_atoms.symbols)):
                    if frame_atoms.symbols[i] == sub[0]:
                        frame_atoms.symbols[i] = sub[1]
        frame_atoms_list.append(frame_atoms)
        frame_atoms.write(f'{file}_povray.{ext}', append=True)

    # For testing purposes
    return frame_atoms_list


def povray_render(atoms, steps, file, ext, atom_subs):

    atom_sub(atoms, atom_subs, steps, file, ext)

    # Allow the user to generate the povray images with reminders of the requirements
    view(atoms)
    print(f'***Crucial Steps***\n'
          f'1. In ASE GUI, navigate to Tools -> Render Scene\n'
          f'2. Change "Output basename" to {file}\n'
          f'3. Select "Render all frames"\n'
          f'4. Deselect "Show output window"\n'
          f'5. Change any other settings (e.g. Atomic texture set) as desired')
    input('***Press Enter to continue once Povray is finished visualising...***')


def gifmaker(steps, file, ext, frames_per_second, pause_time, convert_flags, keep_png_files, test):

    # Generate the list of povray image filenames
    digits = len(str(steps-1))
    indices = [f'%d0{digits}' % i for i in range(steps)]
    filenames = [f'{file}.{index}.png' for index in indices]

    # Get the frames_per_second into a format that ImageMagick accepts for the delay flag
    if frames_per_second <= 1:
        delay = str(1/frames_per_second)
    else:
        delay = f'1x{frames_per_second}'

    # Add in extra first/last frames at the start/end to give a pause of the desired length
    if pause_time is not None:
        pause_frames = int(pause_time * frames_per_second)
        count = 1  # Already one frame in the list
        while count in range(pause_frames):
            filenames.insert(0, f'{file}.{indices[0]}.png')
            filenames.insert(-1, f'{file}.{indices[-1]}.png')
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
    if not test:
        command = (f'convert {convert_options} -delay {delay} %s {file}.gif' % ' '.join(filenames))
        os.system(command)

        # Delete the povray image files if requested
        if not keep_png_files:
            for index in indices:
                os.system(f'rm {file}.{index}.png')
            os.system(f'rm {file}_povray.{ext}')

    # For testing purposes
    return filenames, delay
