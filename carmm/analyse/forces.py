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

    if atoms.calc:
        if 'forces' in atoms.calc.results:
            if not atoms.calc.calculation_required(atoms, ['forces']):
                f = atoms.get_forces()
            if np.amax([np.linalg.norm(f[x]) for x in range(len(atoms))]) <= fmax:
                converged = True
            else:
                converged = False

        if 'stress' in atoms.calc.results:
            if not atoms.calc.calculation_required(atoms, ['stress']):
                s = atoms.get_stress()
            if np.amax([np.linalg.norm(s[x]) for x in range(len(atoms))]) <= fmax:
                converged = True
            else:
                converged = False

    return converged
