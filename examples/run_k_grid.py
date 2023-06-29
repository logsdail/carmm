'''
This example shows how to use the get_k_grid() function which allows to determine reciprocal space sampling
k_grid setting based on simulation cell size and converged value of minimum sampling density.
Useful when running calculations for variable sizes of supercells (e.g. surfaces)
'''

def test_run_k_grid():
    from carmm.run.aims_calculator import get_k_grid
    from ase.build import molecule


    #### Traditional ASE functionality #####
    from data.model_gen import get_example_slab as slab
    slab111 = slab(adsorbate=True)
    slab110 = slab(adsorbate=True, surface="110")
    slab100 = slab(adsorbate=True, surface="100")
    #########
    sampling_density = 0.02 # example value of sampling density /Angstrom
    k_grid111 = get_k_grid(slab111, sampling_density, verbose=True)
    assert k_grid111 == (6, 6, 1)
    k_grid110 = get_k_grid(slab110, sampling_density, verbose=True)
    assert k_grid110 == (5, 6, 1)
    k_grid100 = get_k_grid(slab100, sampling_density, verbose=True)
    assert k_grid100 == (6, 6, 1)
    assert get_k_grid(molecule("CO2"), sampling_density) == None

    # This part tests generating k-grid with the strict definition of reciprocal lattice parameters
    k_grid111 = get_k_grid(slab111, sampling_density, verbose=True, simplify=False)
    assert k_grid111 == (7, 7, 1)
    k_grid110 = get_k_grid(slab110, sampling_density, verbose=True, simplify=False)
    assert k_grid110 == (5, 6, 1)
    k_grid100 = get_k_grid(slab100, sampling_density, verbose=True, simplify=False)
    assert k_grid100 == (6, 6, 1)
    assert get_k_grid(molecule("CO2"), sampling_density, simplify=False) == None

test_run_k_grid()
