import os


def get_aims_calculator(dimensions, spin=None, relativistic=None, k_grid=None, xc="pbe", compute_forces=True, directory='./', **kwargs):
    '''
    Method to return a "default" FHI-aims calculator.
    Note: This file should not be changed without consultation,
          as changes could affect many users in the group.

    Parameters:

        dimensions: Integer
            Determines whether we have a "gas"-phase (0) or "periodic" structure (2 or 3)
        spin: String
            Determines if spin is to be invoked for the calculation
        relativistic: String
            Determines what setting for relativity is to be used
        k_grid: List of integers
            Gives the k-grid sampling in x-, y- and z- direction. e.g. [3, 3, 3]
        xc: String
            XC of choice
        compute_forces: String
            Determines whether forces are enabled ("true") or not enabled ("false").
        directory: String
            ???
        **kwargs:
            Any other keyword arguments that a user wants to set. These are passed through.
        
    Returns:
            FHI_calc: FHI-aims ASE calculator
       
    '''
    from ase.calculators.aims import Aims    

    # Created dictionary to store arguments
    parameter_dict = {}
    
    # Set the XC for the calculation. For LibXC, override_warning_libxc *needs*
    # to be set first, otherwise we get a termination.
    if "libxc" in xc:
        parameter_dict['override_warning_libxc'] = 'true'
    parameter_dict['xc'] = xc

    if dimensions == 2:
        parameter_dict['use_dipole_correction'] = 'true'

    if dimensions >= 2:
        parameter_dict['k_grid'] = k_grid

    if spin is not None:
        parameter_dict['spin'] = 'none'

    if relativistic is None:
        parameter_dict['relativistic'] = ('atomic_zora', 'scalar')

    # Changing to check ASE version, as this determines behaviour of calculator
    from carmm.utils.python_env_check import ase_env_check
    if ase_env_check('3.23.0'):
        # Need a profile for the calculator
        from ase.calculators.aims import AimsProfile
        ase_aims_command = os.environ.get("ASE_AIMS_COMMAND")
        aims_species_dir = os.environ.get("AIMS_SPECIES_DIR")
        if ase_aims_command is None or aims_species_dir is None:
            raise KeyError('Environment variables $ASE_AIMS_COMMAND and $AIMS_SPECIES_DIR are not set')

        fhi_calc = Aims(
            # Load profile from environment variables
            profile=AimsProfile(command=os.environ["ASE_AIMS_COMMAND"],
                                default_species_directory=os.environ["AIMS_SPECIES_DIR"]),
            compute_forces=compute_forces,
            directory=directory,
            # Merged **parameter_dict with **kwargs
            **{**parameter_dict, **kwargs}
        )
    else:
        fhi_calc = Aims(
            compute_forces=compute_forces,
            directory=directory,
            # Merged **parameter_dict with **kwargs
            **{**parameter_dict, **kwargs}
        )

    return fhi_calc


def get_aims_and_sockets_calculator(dimensions,
                                    # i-Pi settings for sockets
                                    port=None, host=None, logfile='socketio.log',
                                    # Debug setting
                                    check_socket=True, verbose=False, codata_warning=True,
                                    # Passthrough of other objects for calculator
                                    **kwargs):
    '''
    Method to return a sockets calculator (for i-Pi based socket connectivity)
    and also an associated FHI-aims calculator for ASE

    Args:
        dimensions: Integer
            Dimensions is _not_ used in this routine, but we make it necessary so
            the passthrough isn't problematic for inexperienced users. See get_aims_calculator()
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
            ASE uses CODATA 2014 and FHI-aims uses CODATA 2002 which yields energy discrepancies.
            The warning message can be turned off if set to False.
    **kwargs:
        Any other keyword arguments that a user wants to set. These are passed through.

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
        host = socket.gethostname()

    # Random port assignment
    if port:
        pass
    else:
        import random
        port = random.randint(12345, 60000)

    if check_socket:
        port = _check_socket(host, port, verbose)

    # **kwargs is a passthrough of keyword arguments
    fhi_calc = get_aims_calculator(dimensions, **kwargs)
    # Add in PIMD command to get sockets working
    fhi_calc.parameters['use_pimd_wrapper']=[host, port]

    # Setup sockets calculator that "wraps" FHI-aims
    from ase.calculators.socketio import SocketIOCalculator
    socket_calc = SocketIOCalculator(fhi_calc, log=logfile, port=port)

    if codata_warning:
        print("You are using i-Pi based socket connectivity between ASE and FHI-aims.")
        print("The communicated energy in Hartree units will be converted to eV in ASE and not FHI-aims.")
        print(
            "The eV/Hartree unit in FHI-aims is given by CODATA 2002 (Web Version 4.0 2003-12-09), Peter J. Mohr, Barry N. Taylor")
        print(
            "ASE uses CODATA 2014, thus the output energy in eV from sockets (ASE conversion) and non-sockets (FHI-aims conversion) will differ.")
        print("e.g. the same atoms.get_total_energy() from ASE will be different using the different calculators.")
        print(
            "Specifically, this is because in non-sockets Aims calculates in a.u., converts to eV (CODATA 2002) and passes that to ASE to output.")
        print(
            "But in sockets, Aims calculates in a.u., passes the a.u. value to ASE through this sockets calculator, ASE converts to eV (CODATA 2014) and outputs that.")
        print(
            "Definition of the constants in each CODATA version can be found at https://wiki.fysik.dtu.dk/ase/_modules/ase/units.html in CODATA{} and create_units().")
        print("PLEASE BE CONSISTENT IN THE UNIT CONVERSION FOR DATA ANALYSIS!")
        print(
            "If you have previous results from the non-sockets calculator, the energy conversion is approximately [Sockets] = [Non-sockets] * 1.00000005204439.")
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
            if verbose: print("Port #" + str(port - 1) + " is unavailable.")
            # Update port
            port += 1
            # Raise issue if port number gets to big!
            if port > 65534:
                raise Exception("No available ports found")
    # Debug statement
    if verbose: print("Port #" + str(port) + " is available.")

    return port


def get_k_grid(model, sampling_density, verbose=False, simple_reciprocal_space_parameters=True):
    """
    Based converged value of reciprocal space sampling provided,
    this function analyses the xyz-dimensions of the simulation cell
    and returns the minimum k-grid as a tuple that can be used
    as a value for the k_grid keyword in the Aims calculator.
    This function allows to maintain consistency in input
    for variable supercell sizes.

    Parameters:
        model: Atoms object
            Periodic model that requires k-grid for calculation in FHI-aims.
        sampling_density: float
            In a simple approach, this is a value of minimum reciprocal space
            sampling required for accuracy of the periodic calculation. The
            value is a fraction between 0 and 1, unit is /Å.
            It is often reported as one k-point per (sampling_density) * 2π Å^-1
            in literature, representing the spacing of k-points along the vector of
            the reciprocal unit cell.

            In FHI-aims, there is a similar parameter called "k_grid_density",
            which is float defining the density of k-point splits along the
            vector of the reciprocal unit cell, unit is nkpts/Å^-1.
            This FHI-aims specific parameter is inversely related to the
            definition above for sampling density as:
            k_grid_density = 1 / (sampling_density * 2 * math.pi)

            Ultimately the choice of definition is a matter of taste, and should be
            explained in any distribution of the outcomes
        simple_reciprocal_space_parameters: bool
            Flag switching between the simplified definition of reciprocal lattice
            parameters, which is 2π/a (a is the real space lattice parameter), and
            the strict definition (see p86 Neil W. Ashcroft Solid State Physics (1976) or
            section 2.4 of https://www.physics-in-a-nutshell.com/article/15/the-reciprocal-lattice).
        verbose: bool
            Flag turning print statements on/off

        Returns:
            float containing 3 integers: (kx, ky, kz)
            or
            None if a non-periodic model is presented
    """

    import math
    import numpy as np

    dimensions = sum(model.pbc)

    if dimensions == 0:
        print("You are studying a molecular system. This has no periodicity, and therefore no k-grid is necessary. "
              "Returning null")
        return None

    # These are lattice vectors
    lattice_v = np.array(model.get_cell())
    # These are lattice parameters
    lattice_param = np.array([np.linalg.norm(v) for v in lattice_v])

    # Check if the model is periodic and with vacuum along a certain axis. There could be vacuum if
    # the lattice parameter is 5 angstrom longer than the range of atomic positions along an axis,
    check_vacuum_and_periodic = np.array((lattice_param - np.ptp(model.get_positions(), axis=0)) > 5) & model.pbc
    if sum(check_vacuum_and_periodic):
        print("There could be vacuum in these axes", np.array(['x', 'y', 'z'])[check_vacuum_and_periodic],
              ", but they are also periodic."
              "\nIf you don't want the model to be treated as periodic in these dimensions,",
              "set pbc for these axes to false, or check if k point sampling is actually 1 in these dimensions")

    if simple_reciprocal_space_parameters:
        # Simplified reciprocal lattice parameters
        reciprocal_param = 2 * math.pi / lattice_param
    else:
        # volume of the cell
        volume = np.dot(lattice_v[0], np.cross(lattice_v[1], lattice_v[2]))
        # These are reciprocal lattice vectors.
        # For definition, see section 2.4 of https://www.physics-in-a-nutshell.com/article/15/the-reciprocal-lattice
        reciprocal_v = [np.cross(lattice_v[(i + 1) % 3], lattice_v[(i + 2) % 3]) * 2 * math.pi / volume
                        for i in range(len(lattice_v))]
        # These are reciprocal lattice parameters
        reciprocal_param = np.array([np.linalg.norm(r_v) for r_v in reciprocal_v])

    k_grid_density = 1 / (sampling_density * 2 * math.pi)
    k_grid = k_grid_density * reciprocal_param
    # Convert k_grid to integer
    k_grid = np.array([math.ceil(k) for k in k_grid])
    # Remove k-sampling if direction is not periodic in any dimension
    k_grid[np.invert(model.pbc)] = 1

    if verbose:
        print("Based on lattice xyz dimensions", "x", round(lattice_param[0], 3), "y", round(lattice_param[1], 3),
              "z", round(lattice_param[2], 3))
        print("and", "one k-point per", str(sampling_density), "* 2π Å^-1",
              "sampling density, the k-grid chosen for periodic calculation is",
              str(k_grid) + ".")
        if not simple_reciprocal_space_parameters:
            print("Please note you are using the strict definition of reciprocal lattice vector here. "
                  "This would generate a slightly denser k-grid than using simple reciprocal space parameters in "
                  "cases where a non-orthogonal cell is used as input.")

    return tuple(k_grid)
