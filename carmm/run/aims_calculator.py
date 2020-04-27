def get_aims_calculator(n):
    '''
    Method to return a "default" FHI-aims calculator.

    TODO: Decided how interactive we want this to be?
    Some of these variables should probably be softcoded

    Parameters:

    n: String
        Determines whether we have a "gas"-phase or "periodic" structure
    '''

    from ase.calculators.aims import Aims

    #"gas" for gas-phase reactants and "periodic" for a periodic systems
    if(n=="gas"):
        return Aims(xc='pbe',
           spin='none',
           vdw_correction_hirshfeld="True",
           relativistic=('atomic_zora','scalar'),
           compute_forces="true"
           )
    else:
        if(n=="periodic"):
            return Aims(xc='pbe',
                spin='none',
                k_grid=(3,3,1),
                vdw_correction_hirshfeld="True",
                relativistic=('atomic_zora','scalar'),
                use_dipole_correction='True',
                compute_forces="true",
                )

    #TODO: What happens if the option is neither of these?