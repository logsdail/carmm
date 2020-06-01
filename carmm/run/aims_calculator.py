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

    if dimensions >= 2:
        fhi_calc.set(k_grid=k_grid)

    return fhi_calc

def get_aims_and_sockets_calculator(dimensions, k_grid=None,
                                    # i-Pi settings for sockets
                                    port=12345, host='localhost', logfile='socketio.log',
                                    # Debug setting
                                    verbose=False):
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
        verbose: Boolean
            For testing of the interface when searching for empty ports.

    Returns:
        Socket_calc: Wrapper for ASE calculator
            Used for i-Pi connectvity, and should be assigned to optimisation/dynamics Object
        FHI_calc: FHI-aims ASE calculator
    '''

    port, port_closed = _check_socket(host, port)
    while port_closed:
        # Debug statement
        if verbose: print("Port #"+str(port-1)+" is closed.")
        # Check next port and update
        port, port_closed = _check_socket(host, port)
        # Raise issue if port number gets to big!
        if port > 65534:
            raise Exception("No available ports found")
    # Debug statement
    if verbose: print("Port #" + str(port) + " is open.")

    fhi_calc = get_aims_calculator(dimensions, k_grid)
    # Add in PIMD command to get sockets working
    fhi_calc.set(use_pimd_wrapper = [host, port])

    # Setup sockets calculator that "wraps" FHI-aims
    from ase.calculators.socketio import SocketIOCalculator
    socket_calc = SocketIOCalculator(fhi_calc, log=logfile, port=port)

    return socket_calc, fhi_calc

def _check_socket(host, port):
    '''
    Function to check if a port is open on the target machine.

    Args:
        host: string
            Name of the host machine on which the port is being tested.
        port: integer
            Port value to test on this iteration

    Returns:
        Integer: Port number as received (if port is open) or updated (if port is closed)
        Boolean: True if port is closed, False if port is open.
    '''
    import socket
    from contextlib import closing

    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        # Check if socket is open. Returns False if so, otherwise True and port is incremented
        if sock.connect_ex((host, port)) == 0:
            return port, False
        else:
            return port+1, True