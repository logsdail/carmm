from ase.io import read, write
from ase.visualize import view
import os

def traj_to_gif(filename, frames_per_second=30, pause_time=0.5, atom_subs=None, keep_temp_files=True):

    '''
    A function which takes a .traj file, visualises it in povray with your desired settings and outputs a .gif file.
    When the function gives you the ase viewer, rotate to your desired view, go to Tools -> Render Scene, select the "Render all frames" option and deselect "Show output window".
    MAKE SURE THAT THE OUTPUT BASENAME IS THE SAME AS THE INITIAL FILENAME (e.g. for 'atoms.traj' the output basename is 'atoms')
    The rest of the settings can be changed at your discretion.
    Then press "Render" and wait for it to finish before closing the window and pressing enter to continue.
    The ImageMagick linux suite is required for this function.
    :param filename: (str) Full name/directory of the .traj file to convert
    :param frames_per_second: (float) Speed of switching images (Default is 30 fps)
    :param pause_time: (float) Time to pause on the first and last images (Default is 0.5 seconds)
    :param atom_subs: (list of lists of strings) Pairs of atomic symbols with the first being changed to the second in all
    images for clearer visualisation
    :param keep_temp_files: (boolean) If set to False, will delete the created povray .traj and .png files after use
    :return: A .gif file of the .traj file, visualised in povray
    '''

    #Retrieve the file name and file extension from (potentially) a full directory path
    if filename == None:
        filename = 'atoms.traj'
    file = filename.split('/')[-1]
    file, ext = file.split('.')

    atoms = read(filename+'@:')
    steps = len(atoms)

    povray_render(atoms, steps, file, ext, atom_subs)

    gifmaker(steps, file, ext, frames_per_second, pause_time, keep_temp_files)

    print("Have a nice day!")


def atom_sub(atoms, atom_subs, steps, file, ext):

    #Replace atoms of one element with another for clearer visualisation
    for frame in range(steps):
        frame_atoms = atoms[frame]
        if atom_subs != None:
            for list in atom_subs:
                for i in range(len(frame_atoms.symbols)):
                    if frame_atoms.symbols[i] == list[0]:
                        frame_atoms.symbols[i] = list[1]
        frame_atoms.write(f'{file}_povray.{ext}', append=True)

    return frame_atoms

def povray_render(atoms, steps, file, ext, atom_subs):

    atom_subs(atoms, atom_subs, steps, file, ext)

    #Allow the user to generate the povray images with reminders of the requirements
    view(atoms)
    print(f'***Remember to change output basename to {file}, select "Render all frames" and deselect "Show output window"***')
    input('Press Enter to continue once Povray is finished visualising...')


def gifmaker(steps, file, ext, frames_per_second, pause_time, keep_png_files):

    #Generate the list of povray image filenames
    digits = len(str(steps-1))
    indices = [f'%0{digits}d' % i for i in range(steps)]
    filenames = [f'{file}.{index}.png' for index in indices]

    #Get the frames_per_second into a format that ImageMagick accepts for the delay flag
    if frames_per_second <= 1:
        delay = str(1/frames_per_second)
    elif frames_per_second > 1:
        delay = f'1x{frames_per_second}'

    #Add in extra first/last frames at the start/end to give a pause of the desired length
    if pause_time != None:
        pause_frames = int(pause_time * frames_per_second)
        count = 1 # Already one frame in the list
        while count in range(pause_frames):
            filenames.insert(0, f'{file}.{indices[0]}.png')
            filenames.insert(-1, f'{file}.{indices[-1]}.png')
            count += 1

    #Execute the ImageMagick convert command in the terminal
    command = (f'convert -verbose -delay {delay} %s -loop 0 {file}.gif' % ' '.join(filenames))
    os.system(command)

    #Delete the povray image files if requested
    if not keep_png_files:
        for index in indices:
            os.system(f'rm {file}.{index}.png')
        os.system(f'{file}_povray.{ext}')

    # For testing purposes
    return filenames, delay
