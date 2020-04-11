"""Helper functions for creating supercells."""

import numpy as np

from ase import Atoms


class SupercellError(Exception):
    """Use if construction of supercell fails"""


def get_deviation_from_optimal_cell_shape(cell, target_shape="sc", norm=None):
    r"""
    Calculates the deviation of the given cell metric from the ideal
    cell metric defining a certain shape. Specifically, the function
    evaluates the expression `\Delta = || Q \mathbf{h} -
    \mathbf{h}_{target}||_2`, where `\mathbf{h}` is the input
    metric (*cell*) and `Q` is a normalization factor (*norm*)
    while the target metric `\mathbf{h}_{target}` (via
    *target_shape*) represent simple cubic ('sc') or face-centered
    cubic ('fcc') cell shapes.

    Parameters:

    cell: 2D array of floats
        Metric given as a (3x3 matrix) of the input structure.
    target_shape: str
        Desired supercell shape. Can be 'sc' for simple cubic or
        'fcc' for face-centered cubic.
    norm: float
        Specify the normalization factor. This is useful to avoid
        recomputing the normalization factor when computing the
        deviation for a series of P matrices.

    """

    if target_shape in ["sc", "simple-cubic"]:
        target_metric = np.eye(3)
    elif target_shape in ["fcc", "face-centered cubic"]:
        target_metric = 0.5 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    if not norm:
        norm = (np.linalg.det(cell) / np.linalg.det(target_metric)) ** (
            -1.0 / 3
        )
    return np.linalg.norm(norm * cell - target_metric)


def find_optimal_cell_shape(
    cell,
    target_size,
    target_shape,
    lower_limit=-2,
    upper_limit=2,
    verbose=False,
):
    """Returns the transformation matrix that produces a supercell
    corresponding to *target_size* unit cells with metric *cell* that
    most closely approximates the shape defined by *target_shape*.

    Parameters:

    cell: 2D array of floats
        Metric given as a (3x3 matrix) of the input structure.
    target_size: integer
        Size of desired super cell in number of unit cells.
    target_shape: str
        Desired supercell shape. Can be 'sc' for simple cubic or
        'fcc' for face-centered cubic.
    lower_limit: int
        Lower limit of search range.
    upper_limit: int
        Upper limit of search range.
    verbose: bool
        Set to True to obtain additional information regarding
        construction of transformation matrix.

    """

    # Set up target metric
    if target_shape in ["sc", "simple-cubic"]:
        target_metric = np.eye(3)
    elif target_shape in ["fcc", "face-centered cubic"]:
        target_metric = 0.5 * np.array(
            [[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float
        )
    if verbose:
        print("target metric (h_target):")
        print(target_metric)

    # Normalize cell metric to reduce computation time during looping
    norm = (
        target_size * np.linalg.det(cell) / np.linalg.det(target_metric)
    ) ** (-1.0 / 3)
    norm_cell = norm * cell
    if verbose:
        print("normalization factor (Q): %g" % norm)

    # Approximate initial P matrix
    ideal_P = np.dot(target_metric, np.linalg.inv(norm_cell))
    if verbose:
        print("idealized transformation matrix:")
        print(ideal_P)
    starting_P = np.array(np.around(ideal_P, 0), dtype=int)
    if verbose:
        print("closest integer transformation matrix (P_0):")
        print(starting_P)

    # Prepare run.
    from itertools import product

    best_score = 1e6
    optimal_P = None
    for dP in product(range(lower_limit, upper_limit + 1), repeat=9):
        dP = np.array(dP, dtype=int).reshape(3, 3)
        P = starting_P + dP
        if int(np.around(np.linalg.det(P), 0)) != target_size:
            continue
        score = get_deviation_from_optimal_cell_shape(
            np.dot(P, norm_cell), target_shape=target_shape, norm=1.0
        )
        if score < best_score:
            best_score = score
            optimal_P = P

    if optimal_P is None:
        print("Failed to find a transformation matrix.")
        return None

    # Finalize.
    if verbose:
        print("smallest score (|Q P h_p - h_target|_2): %f" % best_score)
        print("optimal transformation matrix (P_opt):")
        print(optimal_P)
        print("supercell metric:")
        print(np.round(np.dot(optimal_P, cell), 4))
        print(
            "determinant of optimal transformation matrix: %g"
            % np.linalg.det(optimal_P)
        )
    return optimal_P


def make_supercell(prim, P, wrap=True, tol=1e-5):
    r"""Generate a supercell by applying a general transformation (*P*) to
    the input configuration (*prim*).

    The transformation is described by a 3x3 integer matrix
    `\mathbf{P}`. Specifically, the new cell metric
    `\mathbf{h}` is given in terms of the metric of the input
    configuraton `\mathbf{h}_p` by `\mathbf{P h}_p =
    \mathbf{h}`.

    Parameters:

    prim: ASE Atoms object
        Input configuration.
    P: 3x3 integer matrix
        Transformation matrix `\mathbf{P}`.
    wrap: bool
        wrap in the end
    tol: float
        tolerance for wrapping
    """

    supercell_matrix = P
    supercell = clean_matrix(supercell_matrix @ prim.cell)

    # cartesian lattice points
    lattice_points_frac = lattice_points_in_supercell(supercell_matrix)
    lattice_points = np.dot(lattice_points_frac, supercell)

    superatoms = Atoms(cell=supercell, pbc=prim.pbc)

    for lp in lattice_points:
        shifted_atoms = prim.copy()
        shifted_atoms.positions += lp
        superatoms.extend(shifted_atoms)

    # check number of atoms is correct
    n_target = int(np.round(np.linalg.det(supercell_matrix) * len(prim)))
    if n_target != len(superatoms):
        msg = "Number of atoms in supercell: {}, expected: {}".format(
            n_target, len(superatoms)
        )
        raise SupercellError(msg)

    if wrap:
        superatoms.wrap(eps=tol)

    return superatoms


def lattice_points_in_supercell(supercell_matrix):
    """Find all lattice points contained in a supercell.

    Adapted from pymatgen, which is available under MIT license:
    The MIT License (MIT) Copyright (c) 2011-2012 MIT & The Regents of the
    University of California, through Lawrence Berkeley National Laboratory
    """

    diagonals = np.array(
        [
            [0, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [1, 1, 1],
        ]
    )
    d_points = np.dot(diagonals, supercell_matrix)

    mins = np.min(d_points, axis=0)
    maxes = np.max(d_points, axis=0) + 1

    ar = np.arange(mins[0], maxes[0])[:, None] * np.array([1, 0, 0])[None, :]
    br = np.arange(mins[1], maxes[1])[:, None] * np.array([0, 1, 0])[None, :]
    cr = np.arange(mins[2], maxes[2])[:, None] * np.array([0, 0, 1])[None, :]

    all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
    all_points = all_points.reshape((-1, 3))

    frac_points = np.dot(all_points, np.linalg.inv(supercell_matrix))

    tvects = frac_points[
        np.all(frac_points < 1 - 1e-10, axis=1)
        & np.all(frac_points >= -1e-10, axis=1)
    ]
    assert len(tvects) == round(abs(np.linalg.det(supercell_matrix)))
    return tvects


def clean_matrix(matrix, eps=1e-12):
    """ clean from small values"""
    matrix = np.array(matrix)
    for ij in np.ndindex(matrix.shape):
        if abs(matrix[ij]) < eps:
            matrix[ij] = 0
    return matrix
