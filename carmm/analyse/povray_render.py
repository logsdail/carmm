from ase.io import read, write
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
import subprocess


def povray_render(atoms, output='povray', view=False, atom_subs=None,
                  generic_projection_settings=None, povray_settings=None):
    """
    Saves a .png file showing a povray rendered visualisation of an atoms object.

    Parameters:

    atoms: Atoms object
        Structure to be visualised
    output: String
        Output file basename i.e. 'output'.pov, 'output'.ini etc.
    view: Boolean
        Show the output window
    atom_subs: List of Lists of Strings
        Pairs of atomic symbols with the first being changed to the second for clearer visualisation
        Tip: Find a second atom with a similar atomic radius to the first with a more distinctive colour
    generic_projection_settings: Dictionary
        Settings used by PlottingVariables
        (see https://gitlab.com/ase/ase/-/blob/master/ase/io/utils.py PlottingVariables/__init__ for settings options)
    povray_settings: Dictionary
        Settings used by Povray for visualisation
        (see https://gitlab.com/ase/ase/-/blob/master/ase/io/pov.py POVRAY/__init__ for settings options)

    Returns:

    A .png file of the atoms object, visualised in Povray with the desired settings
    """

    # Default visual settings
    if generic_projection_settings is None:
        generic_projection_settings = {}

    if 'rotation' not in generic_projection_settings:
        generic_projection_settings['rotation'] = '0x,0y,0z'
    if 'radii' not in generic_projection_settings:
        generic_projection_settings['radii'] = 1.0

    if 'colors' not in generic_projection_settings:
        generic_projection_settings['colors'] = None
    else:
        for atom in atoms:
            sym = atom.symbol
            if sym not in generic_projection_settings['colors']:
                generic_projection_settings['colors'][sym] = jmol_colors[atomic_numbers[sym]]

    if povray_settings is None:
        povray_settings = {}

    if 'camera_type' not in povray_settings:
        povray_settings['camera_type'] = 'orthographic angle 5'
    if 'camera_dist' not in povray_settings:
        povray_settings['camera_dist'] = 50

    if view:
        povray_settings['display'] = True
    else:
        povray_settings['display'] = False

    # Camera type information
    if povray_settings['camera_type'] == 'orthographic':
        print('For the orthographic camera type, use "orthographic angle X", where X=5 by default.\n'
              'Increasing/decreasing X has the effect of zooming in/out, respectively.')
    if povray_settings['camera_type'] == 'perspective':
        print('Perspective camera type is much less supported. e.g. Unable to zoom in/out.\n'
              'Instead, try the "orthographic angle X" camera type, where X=5 by default.\n'
              'Increasing/decreasing X has the effect of zooming in/out, respectively.')
    if povray_settings['camera_type'] == 'ultra_wide_angle':
        print('Ultra Wide Angle camera type is much less supported. e.g. Unable to zoom in/out.\n'
              'Instead, try the "orthographic angle X" camera type, where X=5 by default.\n'
              'Increasing/decreasing X has the effect of zooming in/out, respectively.')

    if atom_subs is not None:
        atoms = atom_sub(atoms, atom_subs)

    povobj = write(f'{output}.pov', atoms, **generic_projection_settings, povray_settings=povray_settings)
    try:
        povobj.render()
    except (FileNotFoundError, subprocess.CalledProcessError) as err:
        # Give an error message without stopping (unittest will always fail here)
        print(f'{err}: POVRAY failed to render. Do you have POVRAY installed?')

    # For testing purposes
    return generic_projection_settings, povray_settings


def atom_sub(atoms, atom_subs):
    for sub in atom_subs:
        for atom in range(len(atoms.symbols)):
            if atoms.symbols[atom] == sub[0]:
                atoms.symbols[atom] = sub[1]

    return atoms
