#!/usr/bin/env python3

"""
Example script showing how to add a hydrogen to a zeolite cluster,
then add and rotate an ethanol molecule to the desired position.
"""

def test_adsorbate_placer():
    from ase.io import read
    from ase.build import molecule
    import numpy as np
    from carmm.build.adsorbate_placer import RotationBox
    from ase import Atoms

    mth = molecule('CH3CH2OH')
    site = read("data/H-Y_cluster/H-Y_cluster.xyz")

    h_atom = Atoms('H', positions=[(0, 0, 0)])

    h_placed = RotationBox(h_atom, site, 0, 0, 1.0, lps=2)
    h_placed.place_adsorbate()

    mth_placed = RotationBox(mth, h_placed.ads_and_site, 2, -1, 1.5, lps=1)
    mth_placed.place_adsorbate()

    mth_placed.rotate([-45, 0, -45])

    comp_pos1 = np.array([18.78360766, 5.44998722, -0.87248517])
    comp_pos2 = np.array([20.89866577, 4.39477469, -0.93801803])
    error_pos1 = np.linalg.norm(comp_pos1 - mth_placed.atoms_ads.positions[0], axis=-1)
    error_pos2 = np.linalg.norm(comp_pos2 - mth_placed.atoms_ads.positions[2], axis=-1)

    assert np.isclose(error_pos1, 0, rtol=0, atol=1e-06), f"Error = {error_pos1}"
    assert np.isclose(error_pos2, 0, rtol=0, atol=1e-06), f"Error = {error_pos2}"

    from ase.visualize import view
    view(mth_placed.ads_and_site)

test_adsorbate_placer()