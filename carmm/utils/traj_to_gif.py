from ase.io import read
from ase.visualize import view
from ase.io.trajectory import TrajectoryWriter
import os
from PIL import Image
from carmm.utils.povray_render import povray_render, atom_sub


def traj_to_gif(filename, automatic=False, generic_projection_settings=None, povray_settings=None, frames_per_second=30,
                pause_time=0.5, atom_subs=None, gif_options=None, keep_temp_files=True, **kwargs):
    """
    A function which takes a .traj file, visualises it in povray with your desired settings and outputs a .gif file.

    Parameters:

    filename: String
        Full name/directory of the .traj file to convert
    automatic: Boolean
        If True, automatically renders images using given settings
        If False, opens ASE GUI to allow the user to manually render the images (**FOLLOW THE GIVEN INSTRUCTIONS**)
    generic_projection_settings: Dictionary
        Settings used by PlottingVariables for automatic rendering
        (see https://gitlab.com/ase/ase/-/blob/master/ase/io/utils.py PlottingVariables/__init__ for settings options)
    povray_settings: Dictionary
        Settings used by Povray for automatic rendering
        (see https://gitlab.com/ase/ase/-/blob/master/ase/io/pov.py POVRAY/__init__ for settings options)
    frames_per_second: Float
        Speed of switching images (Default is 30 fps)
    pause_time: Float
        Time (in seconds) to pause on the first and last images (Default is 0.5 seconds)
    atom_subs: List of lists of strings
        Pairs of atomic symbols with the first being changed to the second in
        all images for clearer visualisation. An alternate solution to changing atom colours.
        Tip: Find a second atom with a similar atomic radius to the first with a more distinctive colour
    gif_options: Dictionary of strings
        Settings for the Pillow.Image.save() function. For default setting, don't include
        this parameter. Default options are: "save_all=True, optimize=False, loop=0"
        (see https://pillow.readthedocs.io/en/latest/handbook/image-file-formats.html#gif-saving for full list of options)
    keep_temp_files: Boolean
        If False, will delete the created .png, .pov and .ini files after use. If True, will
        keep these files and produce a .traj file if any atoms were substituted

    Returns:

    A .gif file of the .traj file, visualised in Povray with the desired settings
    """

    # Retrieve the file name and file extension from (potentially) a full directory path
    if filename is None:
        filename = 'atoms.traj'
    file = filename.split('/')[-1]
    ext = file.split('.')[-1]
    file = file.split('.')[:-1]
    file = '.'.join(file)

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
                          generic_projection_settings=generic_projection_settings, povray_settings=povray_settings)
    else:
        if atom_subs is not None:
            for frame in range(steps):
                frame_atoms = atoms[frame]
                atoms[frame] = atom_sub(frame_atoms, atom_subs)
            if keep_temp_files:
                writer = TrajectoryWriter(f'{file}_povray.traj', mode='w')
                for frame in range(steps):
                    writer.write(atoms[frame])
        print(f'***Crucial Steps***\n'
              f'1. In ASE GUI, navigate to Tools -> Render Scene\n'
              f'2. Change "Output basename" to {file}\n'
              f'3. Select "Render all frames"\n'
              f'4. Deselect "Show output window"\n'
              f'5. Change any other settings (e.g. Atomic texture set) as desired')
        view(atoms)
        input('***Press Enter to continue once Povray is finished visualising...***\n')

    gifmaker(file, filenames, frames_per_second, pause_time, gif_options, indices, keep_temp_files)

    print("Happy cooking!")

    # For testing purposes
    return file, ext, steps, atoms, filenames


def gifmaker(file, filenames, frames_per_second, pause_time, gif_options, indices, keep_temp_files):

    duration = [(1 / frames_per_second) * 10**3] * len(filenames)  # Duration in milliseconds

    if pause_time is not None:
        duration[0] = pause_time * 10**3
        duration[-1] = pause_time * 10**3

    # Default gif_options
    if gif_options is None:
        gif_options = {}

    if 'save_all' not in gif_options:
        gif_options['save_all'] = True
    if 'optimize' not in gif_options:
        gif_options['optimize'] = False
    if 'loop' not in gif_options:
        gif_options['loop'] = 0

    # Images
    images = []
    for i in range(len(filenames)):
        try:
            im = Image.open(filenames[i])
            images.append(im)
        except FileNotFoundError as err:
            print(f'{err}')
            break

    try:
        images[0].save(f'{file}.gif', append_images=images[1:], duration=duration, save_all=gif_options['save_all'],
                       optimize=gif_options['optimize'], loop=gif_options['loop'])
    except IndexError as err:
        print(f'{err}: No images')

    # Delete the povray image files if requested
    if not keep_temp_files:
        for index in indices:
            os.system(f'rm {file}.{index}.ini')
            os.system(f'rm {file}.{index}.pov')
            os.system(f'rm {file}.{index}.png')

    # For testing purposes
    return filenames, duration, gif_options
