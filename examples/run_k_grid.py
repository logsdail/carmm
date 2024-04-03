'''
This example shows how to use the get_k_grid() function which allows to determine reciprocal space sampling
k_grid setting based on simulation cell size and converged value of minimum sampling density.
Useful when running calculations for variable sizes of supercells (e.g. surfaces)
'''

def test_run_k_grid():
    from carmm.run.aims_calculator import get_k_grid
    from ase.build import molecule
    from ase import Atoms


    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab111 = slab(adsorbate=True)
    #########
    sampling_density = 0.02 # example value of sampling density /Angstrom
    k_grid111 = get_k_grid(slab111, sampling_density, verbose=True)
    assert k_grid111 == (6, 6, 1)
    assert get_k_grid(molecule("CO2"), sampling_density) == None

    # Test with the strict definition of reciprocal lattice parameters
    k_grid111 = get_k_grid(slab111, sampling_density, verbose=True, simple_reciprocal_space_parameters=False)
    assert k_grid111 == (7, 7, 1)

    # Test for a slab with periodicity in all dimensions
    # This is relevant if using a slab model from pymatgen which will be periodic in 3 dimension by default.
    # The function will proceed as normal, but print a message to warn the user.
    slab111.set_pbc(True)  # Set pbc to True to simulate a pymatgen slab.
    k_grid111 = get_k_grid(slab111, sampling_density, verbose=True, simple_reciprocal_space_parameters=False)
    assert k_grid111 == (7, 7, 3)  # If failed to recognize this as a slab, k_grid111 will be (7, 7, 3)

    # Test for one dimension system
    d = 2.9
    L = 10.0
    wire = Atoms('Au',
                 positions=[[0, L / 2, L / 2]],
                 cell=[d, L, L],
                 pbc=[1, 0, 0])
    k_grid_wire = get_k_grid(wire, sampling_density, verbose=True)
    assert k_grid_wire == (18, 1, 1)

test_run_k_grid()
