def counterpoise_calc(complex_struc, a_id, b_id, fhi_calc=None, a_name='A', b_name='B',
                      verbose=False, dry_run=False):
    """
    This function does counterpoise (CP) correction in one go, assuming a binding complex AB.

    CP correction = A_only + B_only - A_plus_ghost - B_plus_ghost
    A_only has A in the geometry of the binding complex with its own basis
    A_plus_ghost has A in the same geometry as in the complex with B replaced by ghost atoms
    This value should be positive by this definition and should be added to the energy change of interest,
    such as adsorption energy.

    Some references:
    Szalewicz, K., & Jeziorski, B. (1998). The Journal of Chemical Physics, 109(3), 1198–1200.
    https://doi.org/10.1063/1.476667

    Parameters:
        complex_struc: ASE Atoms
            This is the Atoms object which stores the optimized structure of the binding complex
        a_id: list of atom symbols or atom indices for species A
        b_id: list of atom symbols or atom indices for species B
            Please use both symbols or both indices for a_id and b_id.
        fhi_calc: ase.calculators.aims.Aims object
            A FHI_aims calculator constructed with ase Aims
        a_name: optional.
        b_name: optional.
            The name of the two species for your counterpoise correction, which has symbol (or index) a_id and b_id.
        verbose: Flag for printing output.
        dry_run: bool
            Flag for test run (CI-test)

    Returns: float. counterpoise correction value for basis set superposition error
    """
    # Checking if a_id and b_id are mapped correctly
    if not (isinstance(a_id, list) and isinstance(b_id, list)):
        raise TypeError('Please supply a_id and b_id as list.')

    id_type = type(a_id[0])

    for some_id in a_id + b_id:
        if not isinstance(some_id, id_type):
            raise RuntimeError("a_id and b_id should be either lists of indices or lists of strings")

    if id_type is str:
        symbol_not_index = True
        if len(a_id + b_id) == len(complex_struc.symbols):
            raise RuntimeError("The number of symbols are not the same as in the complex")
    elif id_type is int:
        symbol_not_index = False
        if len(a_id + b_id) == len(complex_struc):
            raise RuntimeError("The number of indices are not the same as in the complex")

    species_list = [f'{a_name}_only', f'{a_name}_plus_ghost', f'{b_name}_only', f'{b_name}_plus_ghost']

    energies = []
    # This stores bool lists which indicates if an atom is ghost for each species in the order of
    # ['A_only', 'A_plus_ghost', 'B_only', 'B_plus_ghost']
    ghosts_cp = get_which_are_ghosts(complex_struc, a_id, b_id, symbol_not_index)
    structures_cp = get_structures(complex_struc, a_id, b_id, symbol_not_index)
    for index in range(4):
        if 'compute_forces' in fhi_calc.parameters:
            fhi_calc.parameters.pop('compute_forces')
        fhi_calc.outfilename = species_list[index] + '.out'
        structures_cp[index].calc = fhi_calc
        ghost_calculate(calc=structures_cp[index].calc, atoms=structures_cp[index], ghosts=ghosts_cp[index],
                        dry_run=dry_run)
        energy_i = structures_cp[index].get_potential_energy()
        energies.append(energy_i)

    # This could be done easily without an extra ghost_calculate function if we have the same version of ASE on Gitlab,
    # where ghosts can be specified while constructing the calculator.

    # Counterpoise correction for basis set superposition error
    cp_corr = energies[0] + energies[2] - energies[1] - energies[3]

    if verbose:
        print(species_list, '\n', energies, '\n', cp_corr)

    return cp_corr


def get_which_are_ghosts(complex_struc, a_id, b_id, symbol_not_index):
    """
        This function generates bool lists for writing geometry.in files for counterpoise correction in one go,
        assuming a binding complex AB.
        Parameters:
            complex_struc: ASE Atoms
                This is the Atoms object which stores the optimized structure of the binding complex
            a_id: list of atom symbols or atom indices for species A
            b_id: list of atom symbols or atom indices for species B
                Please use both symbols or both indices for a_id and b_id.
            symbol_not_index: bool
                This indicates whether symbols are used.

        Returns: A list of bool list for each species (Ghost atom is true; Real atom is False).
        """

    # This stores bool lists which indicates if an atom is ghost for each species in the order of
    # ['A_only', 'A_plus_ghost', 'B_only', 'B_plus_ghost']
    if symbol_not_index:
        ghosts_cp = [None, [atom.symbol in b_id for atom in complex_struc],
                     None, [atom.symbol in a_id for atom in complex_struc]]
    else:
        ghosts_cp = [None, [atom.index in b_id for atom in complex_struc],
                     None, [atom.index in a_id for atom in complex_struc]]
    return ghosts_cp


def get_structures(complex_struc, a_id, b_id, symbol_not_index):
    """
        This function generates a list of atoms objects for counterpoise correction, assuming a binding complex AB.
        Parameters:
            complex_struc: ASE Atoms
                This is the Atoms object which stores the optimized structure of the binding complex
            a_id: list of atom symbols or atom indices for species A
            b_id: list of atom symbols or atom indices for species B
                Please use both symbols or both indices for a_id and b_id.
            symbol_not_index: bool
                This indicates whether symbols are used.
        """
    from copy import deepcopy

    # This stores the atoms objects used in cp correction in the order of
    # ['A_only', 'A_plus_ghost', 'B_only', 'B_plus_ghost']
    structures_cp = [deepcopy(complex_struc)] * 4
    # Convert a_id and b_id to index
    if symbol_not_index:
        b_id = [atom.index for atom in complex_struc if atom.symbol in b_id]
        a_id = [atom.index for atom in complex_struc if atom.symbol in a_id]
    del structures_cp[0][b_id]
    del structures_cp[2][a_id]

    return structures_cp


def ghost_calculate(calc, atoms=None, properties=['energy'], system_changes=['positions', 'numbers', 'cell', 'pbc',
                                                                             'initial_charges', 'initial_magmoms'],
                    ghosts=None, dry_run=False):
    """
    This is a modified version of ase.calculators.calculator.FileIOCalculator.calculate to make ghost atoms work
    Args:
        calc: fhi_aims calculator constructed by ase
        atoms: ASE atoms object for counterpoise correction
        properties: list of str. properties to be calculated, default is energy. See original function
        system_changes: list of str. See original function.
        ghosts: bool list. Ghost is Ture and Atom is False. The length is the same as atoms.
        dry_run: flag for CI-test.

    """
    from ase.calculators.calculator import Calculator
    import subprocess
    Calculator.calculate(calc, atoms, properties, system_changes)
    calc.write_input(calc.atoms, properties, system_changes, ghosts=ghosts)
    command = calc.command
    if dry_run:
        command = 'ls'
    subprocess.check_call(command, shell=True, cwd=calc.directory)
    calc.read_results()
