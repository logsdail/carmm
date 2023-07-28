"""
Test script for povray_render and atom_sub functions
Usage example shown first, other tests shown after.
"""


def test_povray_render():
    from carmm.analyse.povray_render import povray_render, atom_sub
    from ase.io import read
    import numpy as np
    from ase.data import covalent_radii

    # Main example script
    atoms = read('data/NH3-H3O_traj/nh3-h3o.traj')

    # An alternative to substituting atoms: changing the colour of specific elements
    # Use 'colors': None for default colours
    color_dict = {
        'H': [255, 255, 255],
        'N': [48, 80, 248],
        'O': [255, 13, 13],
    }
    colors = np.array([color_dict[atom.symbol] for atom in atoms])/255  # RGB colours between 0 and 1

    radius_scale = 0.8
    radius_list = []
    for atomic_number in atoms.get_atomic_numbers():
        if atomic_number == 8:
            radius_list.append(radius_scale * 1)
        else:
            radius_list.append(radius_scale * covalent_radii[atomic_number])

    generic_projection_settings = {
        'rotation': '90x,80y,90z',
        'radii': radius_list,
        'colors': colors,
    }

    povray_settings = {
        'camera_type': 'orthographic angle 7',
        'camera_dist': 100,
        'celllinewidth': 0.07,
        'textures': ['jmol'] * len(atoms),
    }

    povray_render(atoms, generic_projection_settings=generic_projection_settings, povray_settings=povray_settings)

    # Tests

    gen_proj_sett, povray_sett = povray_render(atoms,
                                               generic_projection_settings=generic_projection_settings,
                                               povray_settings=povray_settings)

    assert gen_proj_sett == {
        'rotation': '90x,80y,90z',
        'radii': radius_list,
        'colors': colors,
        'display': False
    }

    assert povray_sett == {
        'camera_type': 'orthographic angle 7',
        'camera_dist': 100,
        'celllinewidth': 0.07,
        'textures': ['jmol'] * len(atoms),
    }

    gen_proj_sett, povray_sett = povray_render(atoms)

    # Default generic_projection_settings and povray_settings
    assert gen_proj_sett == {
        'rotation': '0x,0y,0z',
        'radii': 1.0,
        'colors': None,
    }
    assert povray_sett == {
        'camera_type': 'orthographic angle 5',
        'camera_dist': 50,
    }

    # Atom substitution function, useful for easily changing element colours
    subbed_atoms = atom_sub(atoms, atom_subs=[['N', 'C'], ['O', 'N']])
    assert subbed_atoms.symbols[4] == 'C'
    assert subbed_atoms.symbols[0] == 'N'
    assert len(subbed_atoms) == 8
