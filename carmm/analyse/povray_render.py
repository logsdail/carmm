from ase.io import read, write


def povray_render(atoms, output='povray', view=False, atom_subs=None, generic_projection_settings=None, povray_settings=None):
    """
    Saves a .png file showing a povray rendered visualisation of an atoms object.
    :param atoms: (Atoms object) Structure to be visualised
    :param output: (Str) Output file base
    :param view: (Bool) Show the output window
    :param atom_subs: (List of Lists of Strings) Pairs of atomic symbols with the first being changed to the second
    for clearer visualisation
    :param generic_projection_settings: (Dict) Settings for view of the structure
    (see https://wiki.fysik.dtu.dk/ase/_modules/ase/io/pov.html for settings options)
    :param povray_settings: (Dict) Setting for Povray visualisation
    (see https://wiki.fysik.dtu.dk/ase/_modules/ase/io/pov.html for settings options)
    """
    if generic_projection_settings is None:
        generic_projection_settings = {
            'rotation': '0x,0y,0z',
            'radii': 1.0,
            'colors': None,
        }

    if povray_settings is None:
        povray_settings = {
            'celllinewidth': 0.1,
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
