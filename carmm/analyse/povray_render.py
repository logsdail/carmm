from ase.io import read, write


def povray_render(atoms, output='povray', view=False, atom_subs=None, generic_projection_settings=None, povray_settings=None):
    """
    Saves a .png file showing a povray rendered visualisation of an atoms object.
    :param atoms: (Atoms object) Structure to be visualised
    :param output: (Str) Output file basename i.e. 'output'.pov, 'output'.ini etc.
    :param view: (Bool) Show the output window
    :param atom_subs: (List of Lists of Strings) Pairs of atomic symbols with the first being changed to the second
    for clearer visualisation
    :param generic_projection_settings: (Dict) Settings used by PlottingVariables
    (see https://gitlab.com/ase/ase/-/blob/master/ase/io/utils.py PlottingVariables/__init__ for settings options)
    :param povray_settings: (Dict) Settings used by Povray for visualisation
    (see https://gitlab.com/ase/ase/-/blob/master/ase/io/pov.py POVRAY/__init__ for settings options)
    """
    if generic_projection_settings is None:
        generic_projection_settings = {
            'rotation': '0x,0y,0z',
            'radii': 1.0,
            'colors': None,
        }

    if povray_settings is None:
        povray_settings = {
            'camera_type': 'orthographic angle 5',
            'camera_dist': 50,
        }

    if view:
        povray_settings['display'] = True
    else:
        povray_settings['display'] = False

    if atom_subs is not None:
        atoms = atom_sub(atoms, atom_subs)

    write(f'{output}.pov', atoms, **generic_projection_settings, povray_settings=povray_settings).render()


def atom_sub(atoms, atom_subs):
    for sub in atom_subs:
        for atom in range(len(atoms.symbols)):
            if atoms.symbols[atom] == sub[0]:
                atoms.symbols[atom] = sub[1]

    return atoms
