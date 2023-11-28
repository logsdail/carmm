def counterpoise_calc(complex_struc, a_id, b_id, symbol_not_index, dry_run=False, **kwargs):
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
        dry_run: bool
            Flag for test run (CI-test)
        kwargs: dict
            Keyword Args for get_aims_and_sockets_calculator
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

        socket_calc, fhi_calc = get_aims_and_sockets_calculator(ghost=ghost_list)
        fhi_calc.set(**kwargs)
        fhi_calc.parameters.pop('compute_forces')

        # fhi_calc.outfilename = species+'.out'
        with socket_calc as calc:
            binding_complex.calc = calc
            if dry_run:
                binding_complex.calc = EMT()

            energy = binding_complex.get_potential_energy()
            energies.append(energy)
            subprocess.check_call('mv', 'aims.out', species + '.out')

    # Counterpoise correction for basis set superposition error
    cp_corr = energies[0] + energies[2] - energies[1] - energies[3]

    return cp_corr
