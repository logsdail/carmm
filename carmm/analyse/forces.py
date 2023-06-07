def is_converged(atoms, fmax=0.01):
    '''
    This function takes an atoms object and a force convergence criterion and
    returns True if stored forces are below fmax convergence criterion or False if
    forces are above fmax or there is no stored calculator.
    Useful for multiple geometry optimisations and automated restarts from .traj
    files to avoid unnecessary calculations.

    atoms: Atoms object
    fmax: float
        Convergence criterion in eV/Angstrom, usually 0.01
    '''
    import numpy as np

    converged = False

    if atoms.fhi_calc:
        if not atoms.fhi_calc.calculation_required(atoms, ['forces']):
            f = atoms.get_forces()

            '''
            List comprehension for:
            - Retrieving forces from the calculator
            - Taking their vector norm
            - But only for atoms without constraints
            '''

            if np.amax([np.linalg.norm(f[x]) \
                    for x in range(len(atoms))]) <= fmax:
                converged = True

    return converged
