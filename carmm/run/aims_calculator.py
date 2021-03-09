def get_aims_calculator(dimensions, k_grid=None, xc="pbe", compute_forces="true"):
    '''
    Method to return a "default" FHI-aims calculator.
    Note: This file should not be changed without consultation,
          as changes could affect many users in the group.

    Parameters:

    dimensions: Integer
        Determines whether we have a "gas"-phase (0) or "periodic" structure (2 or 3)
    k_grid: List of integers
        Gives the k-grid sampling in x-, y- and z- direction. e.g. [3, 3, 3]
    xc: String
        XC of choice
    compute_forces: String
        Determines whether forces are enabled ("true") or not enabled ("false").


TODO: Reorder inputs so most necessary are first i.e. xc, compute_forces, k_grid (I think?)
    '''

    from ase.calculators.aims import Aims

    # Default is suitable for molecular calculations
    fhi_calc =  Aims(
                     spin='none',
                     relativistic=('atomic_zora','scalar'),
                     compute_forces=compute_forces
                     )

    # Set the XC for the calculation. For LibXC, override_warning_libxc *needs*
    # to be set first, otherwise we get a termination.
    if "libxc" in xc:
        fhi_calc.set(override_warning_libxc="true")
    fhi_calc.set(xc=xc)

    if dimensions == 2:
        fhi_calc.set(use_dipole_correction='true')

    if dimensions >= 2:
        fhi_calc.set(k_grid=k_grid)

    return fhi_calc

def get_aims_and_sockets_calculator(dimensions, k_grid=None, xc="pbe", compute_forces="true",
                                    # i-Pi settings for sockets
                                    port=None, host=None, logfile='socketio.log',
                                    # Debug setting
                                    check_socket=True, verbose=False, codata_warning=True):
    '''
    Method to return a sockets calculator (for i-Pi based socket connectivity)
    and also an associated FHI-aims calculator for ASE

    Args:
        dimensions: Integer
            See get_aims_calculator()
        k_grid: List of integers
            See get_aims_calculator()
        xc: String
            See get_aims_calculator()
	compute_forces: String
	    See get_aims_calculator()
        port: None or Integer
            The port for connection between FHI-aims and ASE with i-Pi sockets.
            This is fairly arbitrary as long as it doesn't clash with local settings.
            If None an integer between 12345 and 60000 will be picked at random.
        host: String
            Name of host computer for ASE. Necessary for calculations where MPI runs on the compute nodes.
            This is now defaulted to None, and then will self-identify the host name if unidentified.
            Long-term, we should consider removing this variable to insulate the users from having to set it.
        logfile: String
            Location of output from sockets communication.
        check_socket: Boolean
            Whether we want the automatically identify and resolve port clashes.
        verbose: Boolean
            For testing of the interface when searching for empty ports.
        codata_warning: Boolean
            Warn the user about Hartree to eV conversion being performed in ASE rather than FHI-aims.
            ASE uses CODATA 2018 and FHI-aims uses CODATA 2002 which yields energy discrepancies.
            The warning message can be turned off if set to False.

    Returns:
        Socket_calc: Wrapper for ASE calculator
            Used for i-Pi connectvity, and should be assigned to optimisation/dynamics Object
        FHI_calc: FHI-aims ASE calculator
    '''

    import socket

    if host is None:
        # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
        # we need to specifically state what the name of the login node is so the two packages can communicate
        # In order to manage this communication in all situations, here we will find and use the hostname *even*
        # if on the same computer (it should work irrespective)
        host=socket.gethostname()

    # Random port assignment
    if port:
        pass
    else:
        import random
        port = random.randint(12345, 60000)

    if check_socket:
        port = _check_socket(host, port, verbose)

    fhi_calc = get_aims_calculator(dimensions, k_grid, xc, compute_forces)
    # Add in PIMD command to get sockets working
    fhi_calc.set(use_pimd_wrapper = [host, port])

    # Setup sockets calculator that "wraps" FHI-aims
    from ase.calculators.socketio import SocketIOCalculator
    socket_calc = SocketIOCalculator(fhi_calc, log=logfile, port=port)

    if codata_warning:
        print("You are using i-Pi based socket connectivity between ASE and FHI-aims.")
        print("The communicated energy in Hartree units will be converted to eV in ASE and not FHI-aims.")
        print("The eV/Hartree unit in FHI-aims is given by CODATA 2002 (Web Version 4.0 2003-12-09), Peter J. Mohr, Barry N. Taylor")
        print("ASE uses CODATA 2018, thus the energy in eV from ASE and the FHI-aims outputs will differ.")
        print("Please be consistent in the unit conversion for data analysis!")
        print("You can turn off this message by setting 'codata_warning' keyword to False.")

    return socket_calc, fhi_calc

def _check_socket(host, port, verbose=False):
    '''
    Function to check if a port is open on the target machine.

    Args:
        host: string
            Name of the host machine on which the port is being tested.
        port: integer
            Starting port value for testing
        verbose: Boolean
            Whether to output the process of the check on all ports.

    Returns:
        Integer: Port number as received (if port is open) or updated (if port is closed)
    '''
    import socket
    from contextlib import closing

    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        # Check if socket is open. If so, it's unsuitable for use so the port number is updated.
        # Repeat until we find a port number that is not in use currently.
        while not sock.connect_ex((host, port)):
            # Debug statement
            if verbose: print("Port #"+str(port-1)+" is unavailable.")
            # Update port
            port += 1
            # Raise issue if port number gets to big!
            if port > 65534:
                raise Exception("No available ports found")
    # Debug statement
    if verbose: print("Port #" + str(port) + " is available.")

    return port
