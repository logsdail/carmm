def test_geodesic_interpolator():
    from carmm.build.neb.geodesic import GeodesicInterpolator
    from ase.io import read, write
    from ase.build import molecule
    from ase.neb import idpp_interpolate
    import numpy as np

    initial = molecule('C2H6')
    final = molecule('C2H6')
    final[2].position[:-1] = final[5].position[:-1]
    final[4].position[:-1] = final[6].position[:-1]
    final[3].position[:-1] = final[7].position[:-1]

    Geodesic = GeodesicInterpolator(initial, final, 6)
    Geodesic.init_path()

    idpp_interpolate(Geodesic.images)

    ref_images = read("./data/C2H6_path/path.traj", index=":")
    rmsd = [ Geodesic.cart_rmsd(atoms1=ref_image, atoms2=test_im) \
                            for ref_image, test_im in zip(ref_images, Geodesic.images)]

    assert np.all(np.array(rmsd) < 0.01), "Geodesic test failed - difference in geometry greater than 0.01 compared to reference."

test_geodesic_interpolator()

