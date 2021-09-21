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
        if "forces" in atoms.calc.results:
            # TODO: this will probably only work with ase.constraints.FixAtoms,
            #       FixBondLength does not set force to 0
            # extraction of constraints
            constraints = []
            # remove from nested list
            for i in [i.index.tolist() for i in atoms._get_constraints()]:
                for j in i:
                    constraints.append(j)
            '''
            List comprehension for:
            - Retrieving forces from the calculator
            - Taking their vector norm
            - But only for atoms without constraints
            '''
            if np.amax([np.linalg.norm(atoms.calc.results["forces"][x]) \
                    for x in range(len(atoms)) if x not in constraints]) <= fmax:

                converged = True

    return converged
