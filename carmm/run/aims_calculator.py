def get_aims_calculator(dimensions):
    '''
    Method to return a "default" FHI-aims calculator.

    TODO: Some of these variables should probably be softcoded e.g. k-grid

    Parameters:

    dimensions: Integer
        Determines whether we have a "gas"-phase (0) or "periodic" structure (2 or 3)
    '''

    from ase.calculators.aims import Aims

    #"gas" for gas-phase reactants and "periodic" for a periodic systems
    if dimensions == 0:
        return Aims(xc='pbe',
                    spin='none',
                    vdw_correction_hirshfeld="True",
                    relativistic=('atomic_zora','scalar'),
                    compute_forces="true"
                    )
    elif dimensions == 2:
        return Aims(xc='pbe',
                    spin='none',
                    k_grid=(3, 3, 1),
                    vdw_correction_hirshfeld="True",
                    relativistic=('atomic_zora','scalar'),
                    use_dipole_correction='True',
                    compute_forces="true",
                    )
    else: # dimensions == 3:
        return Aims(xc='pbe',
                    spin='none',
                    k_grid=(3, 3, 3),
                    vdw_correction_hirshfeld="True",
                    relativistic=('atomic_zora', 'scalar'),
                    compute_forces="true",
                    )