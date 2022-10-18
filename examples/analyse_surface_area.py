def test_surface_area():
    
    # testing the surface area function
    from carmm.analyse.surface_area import surface_area
    from ase.build import fcc111
    import math

    slab =  fcc111('Pd', (5, 5, 5), vacuum=20.0)
    print(slab)
    surface_slab=surface_area(slab)
    print(surface_slab)
    # Assertion test
    assert(surface_slab == 163.8097876575813)
    
test_surface_area()

