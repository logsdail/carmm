import numpy as np
from ase.geometry.cell import cell_to_cellpar


def get_cell_angles_lengths(cell):
    '''
    Returns cell vectors lengths (a,b,c) as well as different
    angles (alpha, beta, gamma, phi, chi, psi) (in radians).
    '''
    cellpar = cell_to_cellpar(cell)
    cellpar[3:] *= np.pi / 180  # convert angles to radians
    parnames = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    values = {n: p for n, p in zip(parnames, cellpar)}

    volume = abs(np.linalg.det(cell))
    for i, param in enumerate(['phi', 'chi', 'psi']):
        ab = np.linalg.norm(
            np.cross(cell[(i + 1) % 3, :], cell[(i + 2) % 3, :]))
        c = np.linalg.norm(cell[i, :])
        s = np.abs(volume / (ab * c))
        if 1 + 1e-6 > s > 1:
            s = 1.
        values[param] = np.arcsin(s)

    return values


class CellBounds:
    '''
    Class for defining as well as checking limits on
    cell vector lengths and angles

    Parameters:

    bounds: dict
        Any of the following keywords can be used, in
        conjunction with a [low, high] list determining
        the lower and upper bounds:

        a, b, c:
           Minimal and maximal lengths (in Angstrom)
           for the 1st, 2nd and 3rd lattice vectors.

        alpha, beta, gamma:
           Minimal and maximal values (in degrees)
           for the angles between the lattice vectors.

        phi, chi, psi:
           Minimal and maximal values (in degrees)
           for the angles between each lattice vector
           and the plane defined by the other two vectors.

    Example:

    >>> from ase.ga.bulk_utilities import CellBounds
    >>> CellBounds(bounds={'phi': [20, 160],
    ...                    'chi': [60, 120],
    ...                    'psi': [20, 160],
    ...                    'a': [2, 20], 'b': [2, 20], 'c': [2, 20]})
    '''

    def __init__(self, bounds={}):
        self.bounds = {'alpha': [0, np.pi], 'beta': [0, np.pi],
                       'gamma': [0, np.pi], 'phi': [0, np.pi],
                       'chi': [0, np.pi], 'psi': [0, np.pi],
                       'a': [0, 1e6], 'b': [0, 1e6], 'c': [0, 1e6]}

        for param, bound in bounds.items():
            if param not in ['a', 'b', 'c']:
                # Convert angle from degree to radians
                bound = [b * np.pi / 180. for b in bound]
            self.bounds[param] = bound

    def is_within_bounds(self, cell):
        values = get_cell_angles_lengths(cell)
        verdict = True
        for param, bound in self.bounds.items():
            if not (bound[0] <= values[param] <= bound[1]):
                verdict = False
        return verdict


def get_rotation_matrix(u, t):
    '''
    Returns the transformation matrix for rotation over an angle t
    along an axis with direction u.
    '''
    ux, uy, uz = u
    cost, sint = np.cos(t), np.sin(t)
    rotmat = np.array([[(ux**2) * (1 - cost) + cost,
                        ux * uy * (1 - cost) - uz * sint,
                        ux * uz * (1 - cost) + uy * sint],
                       [ux * uy * (1 - cost) + uz * sint,
                        (uy**2) * (1 - cost) + cost,
                        uy * uz * (1 - cost) - ux * sint],
                       [ux * uz * (1 - cost) - uy * sint,
                        uy * uz * (1 - cost) + ux * sint,
                        (uz**2) * (1 - cost) + cost]])
    return rotmat


def convert_for_lammps(atoms):
    """
    Convert a parallel piped (forming right hand basis)
    to lower triangular, low-tilt box that LAMMPS will accept.

    This code draws from a previous LAMMPS interface:
    https://svn.fysik.dtu.dk/projects/ase-extra/trunk/
    ase-extra/trunk/ase/calculators/lammpslib.py
    """
    ase_cell = atoms.get_cell()
    cell = np.matrix.transpose(ase_cell)
    # rotate bases into triangular matrix
    tri_mat = np.zeros((3, 3))
    A = cell[:, 0]
    B = cell[:, 1]
    C = cell[:, 2]
    tri_mat[0, 0] = np.linalg.norm(A)
    Ahat = A / np.linalg.norm(A)
    AxBhat = np.cross(A, B) / np.linalg.norm(np.cross(A, B))
    tri_mat[0, 1] = np.dot(B, Ahat)
    tri_mat[1, 1] = np.linalg.norm(np.cross(Ahat, B))
    tri_mat[0, 2] = np.dot(C, Ahat)
    tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
    tri_mat[2, 2] = np.linalg.norm(np.dot(C, AxBhat))

    atoms.set_cell(tri_mat.T, scale_atoms=True)
    atoms.wrap(pbc=True)

    # "flip" the cell if it is too skewed
    newcell = atoms.get_cell()
    while True:
        xx, yy = newcell[0, 0], newcell[1, 1]
        xy, xz, yz = newcell[1, 0], newcell[2, 0], newcell[2, 1]
        cond1 = 2 * abs(xy) > xx
        cond2 = 2 * abs(xz) > xx
        cond3 = 2 * abs(yz) > yy
        if not cond1 and not cond2 and not cond3:
            break
        if cond1:
            newcell[1, 0] += xx * np.round((0.5 * xx - xy) / xx - 0.5)
        if cond2:
            newcell[2, 0] += xx * np.round((0.5 * xx - xz) / xx - 0.5)
        if cond3:
            newcell[2, 1] += yy * np.round((0.5 * yy - yz) / yy - 0.5)
            newcell[2, 0] += xy * np.round((0.5 * yy - yz) / yy - 0.5)

    atoms.set_cell(newcell, scale_atoms=False)
    atoms.wrap(pbc=True)
