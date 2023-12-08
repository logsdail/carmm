def counterpoise_calc(complex_struc, a_id, b_id, symbol_not_index, fhi_calc=None, a_name='A', b_name='B',
                      verbose=False, dry_run=False):
    """
    This function does counterpoise (CP) correction in one go, assuming a binding complex AB.
    Let's say we have A and B in this complex
    A_only has A in the geometry of the binding complex with its own basis
    A_plus_ghost has A in the same geometry as in the complex with B replaced by ghost atoms
    CP correction = A_only + B_only - A_plus_ghost - B_plus_ghost
    This value should be added to the energy change of interest, such as adsorption energy.

    Parameters:
        complex_struc: ASE Atoms
            This is the Atoms object which stores the optimized structure of the binding complex
        a_id: list of atom symbols or atom indices for species A
        b_id: list of atom symbols or atom indices for species B
            Please use both symbols or both indices for a_id and b_id.
        symbol_not_index: bool
            This indicates whether symbols are used.
        fhi_calc: ase.calculators.aims.Aims object
            A FHI_aims calculator constructed with ase Aims
        a_name: The name of the first species for your counterpoise correction, which has symbol (or index) a_id.
        b_name: The name of the second species for your counterpoise correction, which has symbol (or index) b_id.
        verbose: Flag for printing output.
        dry_run: bool
            Flag for test run (CI-test)

    Returns: counterpoise correction value for basis set superposition error
    """

    species_list = [f'{a_name}_only', f'{a_name}_plus_ghost', f'{b_name}_only', f'{b_name}_plus_ghost']

    energies = []
    ghosts_cp = get_ghosts(complex_struc=complex_struc, a_id=a_id, b_id=b_id, symbol_not_index=symbol_not_index)
    structures_cp = get_structures(complex_struc=complex_struc, a_id=a_id, b_id=b_id, symbol_not_index=symbol_not_index)
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


def get_ghosts(complex_struc, a_id, b_id, symbol_not_index):
    """
        This function generates bool lists for writing ghosts input for counterpoise correction in one go,
        assuming a binding complex AB.
        Parameters:
            complex_struc: ASE Atoms
                This is the Atoms object which stores the optimized structure of the binding complex
            a_id: list of atom symbols or atom indices for species A
            b_id: list of atom symbols or atom indices for species B
                Please use both symbols or both indices for a_id and b_id.
            symbol_not_index: bool
                This indicates whether symbols are used.

        Returns: list of ghosts for each species used in counterpoise correction
        """

    species_list = ['A_only', 'A_plus_ghost', 'B_only', 'B_plus_ghost']
    ghosts_cp = []  # This stores indices of ghost atoms for each species
    for species in species_list:
        binding_complex = complex_struc

        if 'only' in species:
            ghost_list = None
        else:
            if 'A' in species:
                ghost_id = b_id
            elif 'B' in species:
                ghost_id = a_id

            if symbol_not_index:
                ghost_list = [atom.symbol in ghost_id for atom in binding_complex]
            else:
                ghost_list = [atom.index in ghost_id for atom in binding_complex]

        ghosts_cp.append(ghost_list)

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

    species_list = ['A_only', 'A_plus_ghost', 'B_only', 'B_plus_ghost']
    structures_cp = []
    for species in species_list:
        binding_complex = deepcopy(complex_struc)
        binding_complex.set_constraint()
        if 'ghost' in species:
            structures_cp.append(binding_complex)
        elif 'only' in species:

            if 'A' in species:
                del_id = b_id
            elif 'B' in species:
                del_id = a_id

            if symbol_not_index:
                del_list = [atom.index for atom in binding_complex if atom.symbol in del_id]
            else:
                del_list = del_id

            del binding_complex[del_list]

            structures_cp.append(binding_complex)

    return structures_cp


all_changes = ['positions', 'numbers', 'cell', 'pbc',
               'initial_charges', 'initial_magmoms']


def ghost_calculate(calc, atoms=None, properties=['energy'],
                    system_changes=all_changes, ghosts=None, dry_run=False):
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
