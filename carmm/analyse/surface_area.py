def surface_area(atom):
    '''
    A function to calculate the surface of a slab.

    Parameters:

    atom: Atoms object
    testing
    '''
    
    import math

    atom_cell = atom.get_cell_lengths_and_angles()
    surface_area = atom_cell[0] * atom_cell[1] * math.sin(math.radians(atom_cell[
                                                                           5]))  # a, b and the angle in between, this is for one surface, should be multiplied by 2 for surface energy calculations.
    return (surface_area)


