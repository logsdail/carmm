#!/usr/bin/env python3

"""
Example script showing how to add a hydrogen to a zeolite cluster,
then add and rotate an ethanol molecule to the desired position.
"""

def test_adsorbate_placer():
    from ase.io import read
    from ase.build import molecule
    import numpy as np
    from carmm.build.adsorbate_placer import place_adsorbate, rotate_and_place_adsorbate
    from ase import Atoms

    molecule = molecule('CH3CH2OH')

    h_atom = Atoms('H', positions=[(0, 0, 0)])

    site = read("data/H-Y_cluster/H-Y_cluster.xyz")

    zeolite, rotated_ads = place_adsorbate(h_atom, site, 0, 0, 1.0)

    ads_and_site, rotated_ads = rotate_and_place_adsorbate(molecule, zeolite, 1.5,
                                                           2, -1, 0,
                                                           rotation=[-45, 0, -45])

    comp_pos1 = np.array([1.96776122e+01,  5.66864736e+00,  1.89702619e+00])
    comp_pos2 = np.array([1.90289830e+01,  5.04410070e+00,  6.82460837e-01])
    error_pos1 = np.linalg.norm(comp_pos1 - rotated_ads.positions[0], axis=-1)
    error_pos2 = np.linalg.norm(comp_pos2 - rotated_ads.positions[1], axis=-1)

    assert np.isclose(error_pos1, 0, rtol=0, atol=1e-06), f"Error = {error_pos1}"
    assert np.isclose(error_pos2, 0, rtol=0, atol=1e-06), f"Error = {error_pos2}"

    from ase.visualize import view
    view(ads_and_site)

test_adsorbate_placer()
