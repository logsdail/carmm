def surface_area(atoms):
    '''
    A function to calculate the surface of a slab (based on cell vectors)

    Parameters:

    atoms: ASE Atoms object holding slab information
    
    Returns: Surface area, in same units as cell vector lengths (^2)
    '''
    
    import math

    cell_parameters = atoms.cell.cellpar()

    # a, b and the angle in between, this is for one surface, should be multiplied by 2 for surface energy calculations.
    surface_area = cell_parameters[0] * cell_parameters[1] * math.sin(math.radians(cell_parameters[5]))  

    return (surface_area)


