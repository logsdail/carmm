def counterpoise_calc(complex_struc, a_id, b_id, fhi_calc=None, a_name=None, b_name=None,
                      verbose=False, dry_run=False):
    """
    This function does counterpoise (CP) correction for basis set super position error (BSSE) in one go,
    assuming a binding complex AB.

    CP correction = A_only + B_only - A_plus_ghost - B_plus_ghost
    A_only has A in the geometry of the binding complex with its own basis
    A_plus_ghost has A in the same geometry as in the complex with B replaced by ghost atoms
    This value should be positive by this definition and should be added to the energy change of interest,
    such as adsorption energy.

    Some references:
    1. Szalewicz, K., & Jeziorski, B. (1998). The Journal of Chemical Physics, 109(3), 1198–1200.
    https://doi.org/10.1063/1.476667

    2. Galano, A., & Alvarez-Idaboy, J. R. (2006). Journal of Computational Chemistry, 27(11), 1203–1210.
    https://doi.org/10.1002/JCC.20438
    (See the second paragraph in the introduction for a good explanation of why BSSE arises)

    Parameters:
        complex_struc: ASE Atoms
            This is the Atoms object which stores the optimized structure of the binding complex
        a_id: list of atom symbols or atom indices for species A
        b_id: list of atom symbols or atom indices for species B
            Please use both symbols or both indices for a_id and b_id.
        fhi_calc: ase.calculators.aims.Aims object
            A FHI_aims calculator constructed with ase Aims
        a_name: string. optional.
        b_name: string. optional.
            The name of the two species for your counterpoise correction, which has symbol (or index) a_id and b_id.
        verbose: Flag for printing output.
        dry_run: bool
            Flag for test run (CI-test)

    Returns: float. counterpoise correction value for basis set superposition error
    """

    from carmm.utils.python_env_check import ase_env_check

    print("Use version 230612 or newer ones, or empty sites won't work with PBC\n")
    # Check if a_id and b_id are mapped correctly and convert symbols to indices
    a_id, b_id = check_and_convert_id(complex_struc, a_id, b_id)

    # Collect info for input files of A_only, A_plus_ghost, B_only, and B_plus_ghost
    # Get bool lists where a ghost atom is True and a real atom is false
    # Delete B(A) for A(B)_only
    ghosts_lists_cp, structures_cp = gather_info_for_write_input(complex_struc, a_id, b_id)

    # a_name and b_name are default as atoms.symbols
    if a_name is None or b_name is None:
        a_name = str(structures_cp[0].symbols)
        b_name = str(structures_cp[2].symbols)
    # Output names
    species_list = [f'{a_name}_only', f'{a_name}_plus_ghost', f'{b_name}_only', f'{b_name}_plus_ghost']

    # Empty sites does not work with forces. Remove compute_forces.
    if 'compute_forces' in fhi_calc.parameters:
        fhi_calc.parameters.pop('compute_forces')
    if 'sc_accuracy_forces' in fhi_calc.parameters:
        if verbose:
            print('Stop calculation as there is a convergence criterion regarding force.', '\n',
                  'Empty sites does not work with forces. Remove and check how it affects your results.')
        raise KeyError('Found sc_accuracy_forces in parameters! FHI-aims can not calculate forces on empty sites.')
    # Create an empty list to store energies for postprocessing.
    energies = []
    for index in range(4):

        if not ase_env_check('3.23.0'):
            fhi_calc.outfilename = species_list[index] + '.out'
            structures_cp[index].calc = fhi_calc
            # Run the calculation. A workaround. Default calculate function doesn't work with ghost atoms.
            calculate_energy_ghost_compatible(calc=structures_cp[index].calc, atoms=structures_cp[index],
                                              ghosts=ghosts_lists_cp[index], dry_run=dry_run)
            # Get the energy from the converged output.
            energy_i = structures_cp[index].get_potential_energy()
            energies.append(energy_i)
        else:
            fhi_calc.template.outputname = species_list[index] + '.out'
            fhi_calc.parameters['ghosts'] = ghosts_lists_cp[index]
            # Scaled positions does not work with empty sites.
            fhi_calc.parameters['scaled'] = False
            structures_cp[index].calc = fhi_calc
            if dry_run:
                structures_cp[index].calc.template.write_input(fhi_calc.profile, fhi_calc.directory,
                                                               structures_cp[index], fhi_calc.parameters,
                                                               fhi_calc.implemented_properties)
                structures_cp[index].calc.results['energy'] = get_energy_dryrun(fhi_calc.directory,
                                                                                fhi_calc.template.outputname)
                structures_cp[index].calc.atoms = structures_cp[index]

            # Get the energy from the converged output.
            energy_i = structures_cp[index].get_potential_energy()
            energies.append(energy_i)

    # Counterpoise correction for basis set superposition error. See docstring for the formula.
    cp_corr = energies[0] + energies[2] - energies[1] - energies[3]

    if verbose:
        print(species_list, '\n', energies, '\n', cp_corr)

    return cp_corr


def check_and_convert_id(complex_struc, a_id, b_id):
    """
    This function checks if a_id and b_id are supplied correctly (lists of indices or lists of strings),
    and convert symbols to indices.
    Args:
        complex_struc: see counterpoise_calc
        a_id:  see counterpoise_calc
        b_id:  see counterpoise_calc

    Returns:
        a_id, b_id. Two lists of indices for A and B, respectively.

    """
    # Checking if a_id and b_id are mapped correctly
    if not (isinstance(a_id, list) and isinstance(b_id, list)):
        raise TypeError('Please supply a_id and b_id as list.')

    id_type = type(a_id[0])

    for some_id in a_id + b_id:
        if not isinstance(some_id, id_type):
            raise RuntimeError("a_id and b_id should be either lists of indices or lists of strings")

    if id_type is str:
        if len(a_id + b_id) != len(set(complex_struc.symbols)):
            raise RuntimeError("The number of symbols are not the same as in the complex")
        # Convert symbols to indices
        a_id = [atom.index for atom in complex_struc if atom.symbol in a_id]
        b_id = [atom.index for atom in complex_struc if atom.symbol in b_id]
    elif id_type is int:
        if len(a_id + b_id) != len(complex_struc):
            raise RuntimeError("The number of indices are not the same as in the complex")

    return a_id, b_id


def gather_info_for_write_input(complex_struc, a_id, b_id):
    """
    This function collects info for writing input files for A_only, A_plus_ghost, B_only, and B_plus_ghost
    The geometry.in files are written with ase.io.aims.write_aims.
    For species with ghost atoms, write_aims needs a bool list for the keyword argument "ghosts",
    where a ghost atom is True and a normal atom is False.
    Parameters:
        complex_struc: The complex. See counterpoise_calc.
        a_id: indices for A
        b_id: indices for B
    Returns:
        A list of four bool lists (Ghost atom is True; Real atom is False)
        A list of atoms objects
        Both in the order of A_only, A_plus_ghost, B_only, B_plus_ghost.
    """

    # Determine whether an atom is ghost atom or not.
    ghosts_cp = [None, [atom.index in b_id for atom in complex_struc],
                 None, [atom.index in a_id for atom in complex_struc]]
    # A list of four bool lists. For ?_only, the value is None.
    # For ?_plus_ghost, the bool list has the same length as complex_struc. (Ghost atom is True; Real atom is False)
    # Order: A_only, A_plus_ghost, B_only, B_plus_ghost

    from copy import deepcopy
    # Prepare atoms objects. Delete A or B as appropriate
    a_only = deepcopy(complex_struc)
    del a_only[b_id]  # Delete B from the complex.
    b_only = deepcopy(complex_struc)
    del b_only[a_id]  # Delete A from the complex.
    structures_cp = [a_only, complex_struc, b_only, complex_struc]
    # Order: A_only, A_plus_ghost, B_only, B_plus_ghost
    # ?_plus_ghost use the same atoms object as the complex

    return ghosts_cp, structures_cp


def calculate_energy_ghost_compatible(calc, atoms=None, properties=['energy'],
                                      system_changes=['positions', 'numbers', 'cell', 'pbc',
                                                      'initial_charges', 'initial_magmoms'],
                                      ghosts=None, dry_run=False):
    """
    This is a modified version of ase.calculators.calculator.FileIOCalculator.calculate to make ghost atoms work
    Do the calculation and read the results.

    This is a workaround. The same could be done easily if we were using the same version of ASE on Gitlab.
    The Aims calculator were rewritten where ghosts can be specified while constructing the calculator
    (not available in current release)

    Args:
        calc: fhi_aims calculator constructed by ASE
        atoms: ASE atoms object
        properties: list of str. properties to be calculated, default is energy. See original function
        system_changes: list of str. See original function.
        ghosts: bool list. Ghost is Ture and Atom is False. The length is the same as atoms.
        dry_run: flag for CI-test.

    """
    from ase.calculators.calculator import Calculator
    import subprocess, os
    Calculator.calculate(calc, atoms, properties, system_changes)
    # Write inputfiles. Scaled positions does not work with empty sites.
    calc.write_input(calc.atoms, properties, system_changes, ghosts=ghosts, scaled=False)
    command = calc.command

    if dry_run:  # Only for CI tests
        command = ''  # Used to be 'ls'
    converged = False
    if os.path.exists(calc.directory+'/'+calc.outfilename):
        converged = calc.read_convergence()
    if (not converged) or dry_run:
        subprocess.check_call(command, shell=True, cwd=calc.directory)

    calc.read_results()


# Lazy work around
def get_energy_dryrun(dir, outputname):
    """Parse the energy from the aims.out file"""
    import numpy as np
    f = open(dir / outputname, 'r')
    fd = f.readlines()
    for line in fd:
        if 'Total energy corrected' in line:
            energy_line = line

    return float(energy_line.split()[5])
