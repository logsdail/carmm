def get_aims_calculator(dimensions, k_grid=None):
    '''
    Method to return a "default" FHI-aims calculator.
    Note: This file should not be changed without consultation,
          as changes could affect many users in the group.

    TODO: Some of these variables should probably be softcoded e.g. k-grid

    Parameters:

    dimensions: Integer
        Determines whether we have a "gas"-phase (0) or "periodic" structure (2 or 3)
    k_grid: List of integers
        Gives the k-grid sampling in x-, y- and z- direction. e.g. [3, 3, 3]
    '''

    from ase.calculators.aims import Aims

    # Default is suitable for molecular calculations
    fhi_calc =  Aims(xc='pbe',
                     spin='none',
                     relativistic=('atomic_zora','scalar'),
                     compute_forces="true"
                     )

    if dimensions == 2:
        fhi_calc.set(use_dipole_correction='True')
        fhi_calc.set(k_grid=k_grid)

    if dimensions > 2:
        fhi_calc.set(k_grid=k_grid)

    return fhi_calc

def get_aims_and_sockets_calculator(dimensions, k_grid=None, port=12345, host='localhost', logfile='socketio.log'):
    '''
    Method to return a sockets calculator (for i-Pi based socket connectivity)
    and also an associated FHI-aims calculator for ASE

    Args:
        dimensions: Integer
            See get_aims_calculator()
        k_grid: List of integers
            See get_aims_calculator()
        port: Integer
            The port for connection between FHI-aims and ASE with i-Pi sockets.
            This is fairly arbitrary as long as it doesn't clash with local settings.
        host: String
            Name of host computer for ASE. Necessary for calculations where MPI runs on the compute nodes

    Returns:
        Socket_calc: Wrapper for ASE calculator
            Used for i-Pi connectvity, and should be assigned to optimisation/dynamics Object
        FHI_calc: FHI-aims ASE calculator
    '''

    fhi_calc = get_aims_calculator(dimensions, k_grid)
    # Add in PIMD command to get sockets working
    fhi_calc.set(use_pimd_wrapper = [host, port])

    # Setup sockets calculator that "wraps" FHI-aims
    from ase.calculators.socketio import SocketIOCalculator
    socket_calc = SocketIOCalculator(fhi_calc, log=logfile, port=port)

    return socket_calc, fhi_calc
