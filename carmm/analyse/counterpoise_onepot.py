def counterpoise_calc(complex_struc, a_id, b_id, symbol_not_index, fhi_calc=None, dry_run=False):
    """
    This function does counterpoise correction in one go, assuming a binding complex AB.
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
        dry_run: bool
            Flag for test run (CI-test)

    Returns: counterpoise correction value for basis set superposition error
    """
    import subprocess
    from copy import deepcopy
    from carmm.run.aims_calculator import get_aims_and_sockets_calculator
    from ase.calculators.emt import EMT
    # Let's say we have A and B in this complex
    # ?_only has A or B in the geometry of the binding complex with its own basis
    # ?_plus_ghost has A or B in the same geometry with ghost atoms added
    species_list = ['A_only', 'A_plus_ghost', 'B_only', 'B_plus_ghost']
    print(species_list)
    energies = []
    for species in species_list:
        binding_complex = deepcopy(complex_struc)
        binding_complex.set_constraint()
        if 'A' in species:
            ghost_id = b_id
        elif 'B' in species:
            ghost_id = a_id

        if symbol_not_index:
            del_list = [atom.index for atom in binding_complex if atom.symbol in ghost_id]
            ghost_list = [atom.symbol in ghost_id for atom in binding_complex]
        else:
            del_list = ghost_id
            ghost_list = [atom.index in ghost_id for atom in binding_complex]

        if 'only' in species:
            del binding_complex[del_list]
            ghost_list = None
        """
        This could be done easily without an extra ghost_calculate function if we have the same version of ASE on Gitlab
        socket_calc, fhi_calc = get_aims_and_sockets_calculator(ghosts=ghost_list, dimensions=dimensions)
        fhi_calc.set(**kwargs)
        fhi_calc.parameters.pop('compute_forces')
        fhi_calc.outputname = species + '.out'
        with socket_calc as calc:
            binding_complex.calc = calc
            if dry_run:
                binding_complex.calc = EMT()

            energy = binding_complex.get_potential_energy()
            energies.append(energy)
        """
        if 'compute_forces' in fhi_calc.parameters:
            fhi_calc.parameters.pop('compute_forces')
        fhi_calc.outfilename = species + '.out'
        binding_complex.calc = fhi_calc
        ghost_calculate(calc=binding_complex.calc, atoms=binding_complex, ghosts=ghost_list, dry_run=dry_run)
        energy = binding_complex.get_potential_energy()
        energies.append(energy)
        print(energies)
    # Counterpoise correction for basis set superposition error
    cp_corr = energies[0] + energies[2] - energies[1] - energies[3]
    print(energies)
    return cp_corr


all_changes = ['positions', 'numbers', 'cell', 'pbc',
               'initial_charges', 'initial_magmoms']


def ghost_calculate(calc, atoms=None, properties=['energy'],
                    system_changes=all_changes, ghosts=None, dry_run=False):
    """
    This is a modified version of ase.calculators.calculator.FileIOCalculator.calculate to make ghost atoms work
    Args:
        calc:
        atoms:
        properties:
        system_changes:
        ghosts:
        dry_run:

    """
    from ase.calculators.calculator import Calculator
    from ase.calculators.emt import EMT
    import subprocess
    Calculator.calculate(calc, atoms, properties, system_changes)
    calc.write_input(calc.atoms, properties, system_changes, ghosts=ghosts)
    command = calc.command
    if dry_run:
        command = 'ls'
    subprocess.check_call(command, shell=True, cwd=calc.directory)
    calc.read_results()
