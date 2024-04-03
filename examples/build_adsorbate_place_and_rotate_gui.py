#!/usr/bin/env python3

"""
Example script showing how to add a hydrogen to a zeolite cluster,
then add and rotate an ethanol molecule to the desired position.
"""

def test_adsorbate_placer_gui():
    from ase.io import read
    from ase.build import molecule
    import numpy as np
    from carmm.build.adsorbate_placer import RotationBox
    from carmm.build.adsorbate_placer_gui import rotationbox_gui
    from ase import Atoms

    for cutoff_mult in [1, 1.2]:

        mth = molecule('CH3CH2OH')
        site = read("data/H-Y_cluster/H-Y_cluster.xyz")

        h_atom = Atoms('H', positions=[(0, 0, 0)])

        h_placed = RotationBox(h_atom, site, 0, 0, 1.0, lps=2)
        h_placed.place_adsorbate()

        mth_placed = RotationBox(mth, h_placed.ads_and_site, 2, -1, 1.5, lps=1, cutoff_mult=cutoff_mult)
        mth_placed.place_adsorbate()

    try:
        rotationbox_gui(mth_placed)
    except Exception as error:
        print(f"Error in rotation box GUI test, {type(error).__name__}")

test_adsorbate_placer_gui()
