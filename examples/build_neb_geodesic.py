def test_geodesic_interpolator():
    from carmm.build.neb.geodesic import GeodesicInterpolator
    from ase.io import read
    import numpy as np

    initial = read("./data/H+CH4_CH3+H2_path/H+CH4_CH3+H2.xyz", index=0)
    final = read("./data/H+CH4_CH3+H2_path/H+CH4_CH3+H2.xyz", index=-1)

    # Test the initialised path
    Geodesic = GeodesicInterpolator(initial, final, 6)
    Geodesic.init_path()
    ref_images_initial = read("./data/H+CH4_CH3+H2_path/ref_path_initial.traj", index=":")
    rmsd_initial = [ Geodesic.cart_rmsd(atoms1=ref_image, atoms2=test_im) \
                            for ref_image, test_im in zip(ref_images_initial, Geodesic.images)]
    assert np.all(np.array(rmsd_initial) < 0.01), \
        "Geodesic test failed - difference in geometry greater than 0.01 compared to reference."

    # Test the iteratively smoothed path
    Geodesic.sweep_iterative(sweeperiter=5)
    ref_images_sweep = read("./data/H+CH4_CH3+H2_path/ref_path_sweep.traj", index=":")
    rmsd_sweep = [ Geodesic.cart_rmsd(atoms1=ref_image, atoms2=test_im) \
                            for ref_image, test_im in zip(ref_images_sweep, Geodesic.images)]
    assert np.all(np.array(rmsd_sweep) < 0.01), \
        "Geodesic test failed - difference in geometry greater than 0.01 compared to reference."

test_geodesic_interpolator()

