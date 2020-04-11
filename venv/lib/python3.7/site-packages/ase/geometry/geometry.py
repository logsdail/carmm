# Copyright (C) 2010, Jesper Friis
# (see accompanying license files for details).

"""Utility tools for atoms/geometry manipulations.
   - convenient creation of slabs and interfaces of
different orientations.
   - detection of duplicate atoms / atoms within cutoff radius
"""

import itertools
import numpy as np
from ase.geometry import complete_cell
from ase.geometry.minkowski_reduction import minkowski_reduce
from ase.utils import pbc2pbc


def translate_pretty(fractional, pbc):
    """Translates atoms such that fractional positions are minimized."""

    for i in range(3):
        if not pbc[i]:
            continue

        indices = np.argsort(fractional[:, i])
        sp = fractional[indices, i]

        widths = (np.roll(sp, 1) - sp) % 1.0
        fractional[:, i] -= sp[np.argmin(widths)]
        fractional[:, i] %= 1.0
    return fractional


def wrap_positions(positions, cell, pbc=True, center=(0.5, 0.5, 0.5),
                   pretty_translation=False, eps=1e-7):
    """Wrap positions to unit cell.

    Returns positions changed by a multiple of the unit cell vectors to
    fit inside the space spanned by these vectors.  See also the
    :meth:`ase.Atoms.wrap` method.

    Parameters:

    positions: float ndarray of shape (n, 3)
        Positions of the atoms
    cell: float ndarray of shape (3, 3)
        Unit cell vectors.
    pbc: one or 3 bool
        For each axis in the unit cell decides whether the positions
        will be moved along this axis.
    center: three float
        The positons in fractional coordinates that the new positions
        will be nearest possible to.
    pretty_translation: bool
        Translates atoms such that fractional coordinates are minimized.
    eps: float
        Small number to prevent slightly negative coordinates from being
        wrapped.

    Example:

    >>> from ase.geometry import wrap_positions
    >>> wrap_positions([[-0.1, 1.01, -0.5]],
    ...                [[1, 0, 0], [0, 1, 0], [0, 0, 4]],
    ...                pbc=[1, 1, 0])
    array([[ 0.9 ,  0.01, -0.5 ]])
    """

    if not hasattr(center, '__len__'):
        center = (center,) * 3

    pbc = pbc2pbc(pbc)
    shift = np.asarray(center) - 0.5 - eps

    # Don't change coordinates when pbc is False
    shift[np.logical_not(pbc)] = 0.0

    assert np.asarray(cell)[np.asarray(pbc)].any(axis=1).all(), (cell, pbc)

    cell = complete_cell(cell)
    fractional = np.linalg.solve(cell.T,
                                 np.asarray(positions).T).T - shift

    if pretty_translation:
        fractional = translate_pretty(fractional, pbc)
        shift = np.asarray(center) - 0.5
        shift[np.logical_not(pbc)] = 0.0
        fractional += shift
    else:
        for i, periodic in enumerate(pbc):
            if periodic:
                fractional[:, i] %= 1.0
                fractional[:, i] += shift[i]

    return np.dot(fractional, cell)


def get_layers(atoms, miller, tolerance=0.001):
    """Returns two arrays describing which layer each atom belongs
    to and the distance between the layers and origo.

    Parameters:

    miller: 3 integers
        The Miller indices of the planes. Actually, any direction
        in reciprocal space works, so if a and b are two float
        vectors spanning an atomic plane, you can get all layers
        parallel to this with miller=np.cross(a,b).
    tolerance: float
        The maximum distance in Angstrom along the plane normal for
        counting two atoms as belonging to the same plane.

    Returns:

    tags: array of integres
        Array of layer indices for each atom.
    levels: array of floats
        Array of distances in Angstrom from each layer to origo.

    Example:

    >>> import numpy as np
    >>> from ase.spacegroup import crystal
    >>> atoms = crystal('Al', [(0,0,0)], spacegroup=225, cellpar=4.05)
    >>> np.round(atoms.positions, decimals=5)
    array([[ 0.   ,  0.   ,  0.   ],
           [ 0.   ,  2.025,  2.025],
           [ 2.025,  0.   ,  2.025],
           [ 2.025,  2.025,  0.   ]])
    >>> get_layers(atoms, (0,0,1))  # doctest: +ELLIPSIS
    (array([0, 1, 1, 0]...), array([ 0.   ,  2.025]))
    """
    miller = np.asarray(miller)

    metric = np.dot(atoms.cell, atoms.cell.T)
    c = np.linalg.solve(metric.T, miller.T).T
    miller_norm = np.sqrt(np.dot(c, miller))
    d = np.dot(atoms.get_scaled_positions(), miller) / miller_norm

    keys = np.argsort(d)
    ikeys = np.argsort(keys)
    mask = np.concatenate(([True], np.diff(d[keys]) > tolerance))
    tags = np.cumsum(mask)[ikeys]
    if tags.min() == 1:
        tags -= 1

    levels = d[keys][mask]
    return tags, levels


def find_mic(v, cell, pbc=True):
    """Finds the minimum-image representation of vector(s) v"""

    pbc = cell.any(1) & pbc2pbc(pbc)
    v = np.array(v)
    single = len(v.shape) == 1
    v = np.atleast_2d(v)

    if np.sum(pbc) > 0:
        cell = complete_cell(cell)
        rcell, _ = minkowski_reduce(cell, pbc=pbc)

        # in a Minkowski-reduced cell we only need to test nearest neighbors
        cs = [np.arange(-1 * p, p + 1) for p in pbc]
        neighbor_cells = list(itertools.product(*cs))

        positions = wrap_positions(v, rcell, pbc=pbc, eps=0)
        vmin = positions.copy()
        vlen = np.linalg.norm(positions, axis=1)
        for nbr in neighbor_cells:
            trial = positions + np.dot(rcell.T, nbr)
            trial_len = np.linalg.norm(trial, axis=1)

            indices = np.where(trial_len < vlen)
            vmin[indices] = trial[indices]
            vlen[indices] = trial_len[indices]
    else:
        vmin = v.copy()
        vlen = np.linalg.norm(vmin, axis=1)

    if single:
        return vmin[0], vlen[0]
    else:
        return vmin, vlen


def get_angles(v1, v2, cell=None, pbc=None):
    """Get angles formed by two lists of vectors.

    calculate angle in degrees between vectors v1 and v2

    Set a cell and pbc to enable minimum image
    convention, otherwise angles are taken as-is.
    """

    # Check if using mic
    if cell is not None or pbc is not None:
        if cell is None or pbc is None:
            raise ValueError("cell or pbc must be both set or both be None")

        v1 = find_mic(v1, cell, pbc)[0]
        v2 = find_mic(v2, cell, pbc)[0]

    nv1 = np.linalg.norm(v1, axis=1)[:, np.newaxis]
    nv2 = np.linalg.norm(v2, axis=1)[:, np.newaxis]
    if (nv1 <= 0).any() or (nv2 <= 0).any():
        raise ZeroDivisionError('Undefined angle')
    v1 /= nv1
    v2 /= nv2

    # We just normalized the vectors, but in some cases we can get
    # bad things like 1+2e-16.  These we clip away:
    angles = np.arccos(np.einsum('ij,ij->i', v1, v2).clip(-1.0, 1.0))
    return np.degrees(angles)


def get_distances(p1, p2=None, cell=None, pbc=None):
    """Return distance matrix of every position in p1 with every position in p2

    if p2 is not set, it is assumed that distances between all positions in p1
    are desired. p2 will be set to p1 in this case.

    Use set cell and pbc to use the minimum image convention.
    """
    p1 = np.atleast_2d(p1)
    if p2 is None:
        np1 = len(p1)
        ind1, ind2 = np.triu_indices(np1, k=1)
        D = p1[ind2] - p1[ind1]
    else:
        p2 = np.atleast_2d(p2)
        D = (p2[np.newaxis, :, :] - p1[:, np.newaxis, :]).reshape((-1, 3))

    # Check if using mic
    if cell is not None or pbc is not None:
        if cell is None or pbc is None:
            raise ValueError("cell or pbc must be both set or both be None")

        D, D_len = find_mic(D, cell, pbc)
    else:
        D_len = np.sqrt((D**2).sum(1))

    if p2 is None:
        Dout = np.zeros((np1, np1, 3))
        Dout[(ind1, ind2)] = D
        Dout -= np.transpose(Dout, axes=(1, 0, 2))

        Dout_len = np.zeros((np1, np1))
        Dout_len[(ind1, ind2)] = D_len
        Dout_len += Dout_len.T
        return Dout, Dout_len

    # Expand back to matrix indexing
    D.shape = (-1, len(p2), 3)
    D_len.shape = (-1, len(p2))

    return D, D_len


def get_duplicate_atoms(atoms, cutoff=0.1, delete=False):
    """Get list of duplicate atoms and delete them if requested.

    Identify all atoms which lie within the cutoff radius of each other.
    Delete one set of them if delete == True.
    """
    from scipy.spatial.distance import pdist
    dists = pdist(atoms.get_positions(), 'sqeuclidean')
    dup = np.nonzero(dists < cutoff**2)
    rem = np.array(_row_col_from_pdist(len(atoms), dup[0]))
    if delete:
        if rem.size != 0:
            del atoms[rem[:, 0]]
    else:
        return rem


def _row_col_from_pdist(dim, i):
    """Calculate the i,j index in the square matrix for an index in a
    condensed (triangular) matrix.
    """
    i = np.array(i)
    b = 1 - 2 * dim
    x = (np.floor((-b - np.sqrt(b**2 - 8 * i)) / 2)).astype(int)
    y = (i + x * (b + x + 2) / 2 + 1).astype(int)
    if i.shape:
        return list(zip(x, y))
    else:
        return [(x, y)]


def permute_axes(atoms, permutation):
    """Permute axes of unit cell and atom positions. Considers only cell and
    atomic positions. Other vector quantities such as momenta are not
    modified."""
    assert (np.sort(permutation) == np.arange(3)).all()

    permuted = atoms.copy()
    scaled = permuted.get_scaled_positions()
    permuted.set_cell(permuted.cell.permute_axes(permutation),
                      scale_atoms=False)
    permuted.set_scaled_positions(scaled[:, permutation])
    permuted.set_pbc(permuted.pbc[permutation])
    return permuted
