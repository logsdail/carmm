def is_converged(atoms, fmax):
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
        if np.amax([np.linalg.norm(f) for f in atoms.get_forces()]) <= fmax:
            converged = True

    return converged
