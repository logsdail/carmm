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

    ads_and_site, rotated_ads = rotate_and_place_adsorbate(molecule, zeolite, 1.0,
                                                           2, 0, 0,
                                                           rotation=[45, 0, -45])

    comp_pos1 = np.array([2.00311656e+01, 5.31509397e+00, 1.89702619e+00])
    comp_pos2 = np.array([2.04405930e+01, 4.53669268e+00, 2.54385921e+00])

    error_pos1 = np.linalg.norm(comp_pos1 - rotated_ads.positions[0], axis=-1)
    error_pos2 = np.linalg.norm(comp_pos2 - rotated_ads.positions[-1], axis=-1)

    assert (np.isclose(error_pos1, 0, rtol=0, atol=1e-06))
    assert (np.isclose(error_pos2, 0, rtol=0, atol=1e-06))

test_adsorbate_placer()
