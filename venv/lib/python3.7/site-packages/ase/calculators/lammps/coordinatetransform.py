import numpy as np
from ase.geometry import wrap_positions

# Order in which off-diagonal elements are checked for strong tilt
FLIP_ORDER = [(1, 0, 0), (2, 0, 0), (2, 1, 1)]


class Prism:
    """The representation of the unit cell in LAMMPS

    The main purpose of the prism-object is to create suitable
    string representations of prism limits and atom positions
    within the prism.

    :param cell: cell in ase coordinate system
    :param pbc: periodic boundaries
    :param tolerance: precision for skewness test

    **Implementation**

    LAMMPS requires triangular matrixes without a strong tilt.
    Therefore the 'Prism'-object contains three coordinate systems:

    - ase_cell (the simulated system in the ASE coordination system)
    - lammps_tilt (ase-cell rotated to be an lower triangular matrix)
    - lammps_cell (same volume as tilted cell, but reduce edge length)

    The translation between 'ase_cell' and 'lammps_cell' is down with a
    rotation matrix 'rot_mat' obtained from a QR decomposition.  The
    transformation between 'lammps_tilt' and 'lammps_cell' is done by
    changing the off-diagonal elements.  'flip' saves the modified
    elements for later reversal.  All vectors except positions are just
    rotated with 'rot_mat'.  Positions are rotated and wrapped into
    'lammps_cell'.  Translating results back from LAMMPS needs first the
    unwrapping of positions.  Then all vectors and the unit-cell are
    rotated back into the ASE coordination system.  This can fail as
    depending on the simulation run LAMMPS might have changed the
    simulation box significantly.  This is for example a problem with
    hexagonal cells.  LAMMPS might also wrap atoms across periodic
    boundaries, which can lead to problems for example NEB
    calculations.
    """

    # !TODO: derive tolerence from cell-dimensions
    def __init__(self, cell, pbc=(True, True, True), tolerance=1.0e-8):
        # Use QR decomposition to get the lammps cell
        #    rot_mat * lammps_tilt^T = ase_cell^T
        # => lammps_tilt * rot_mat^T = ase_cell
        # => lammps_tilt             = ase_cell * rot_mat
        qtrans, ltrans = np.linalg.qr(cell.T, mode="complete")
        self.rot_mat = qtrans
        self.lammps_tilt = ltrans.T
        self.ase_cell = cell
        self.tolerance = tolerance
        self.pbc = pbc

        # LAMMPS requires positive values on the diagonal of the
        # triangluar matrix -> mirror if necessary
        for i in range(3):
            if self.lammps_tilt[i][i] < 0.0:
                mirror_mat = np.eye(3)
                mirror_mat[i][i] = -1.0
                self.lammps_tilt = np.dot(mirror_mat, self.lammps_tilt.T).T
                self.rot_mat = np.dot(self.rot_mat, mirror_mat)

        self.lammps_cell = self.lammps_tilt.copy()

        # LAMMPS minimizes the edge length of the parallelepiped
        # What is ment with 'flip': cell 2 is transformed into cell 1
        # cell 2 = 'lammps_tilt'; cell 1 = 'lammps_cell'
        # o-----------------------------/==o-----------------------------/--o
        #  \                        /--/    \                        /--/
        #   \                   /--/         \                   /--/
        #    \         1    /--/              \   2          /--/
        #     \         /--/                   \         /--/
        #      \    /--/                        \    /--/
        #       o==/-----------------------------o--/
        # !TODO: handle extreme tilt (= off-diagonal > 1.5)
        self.flip = np.array(
            [
                abs(self.lammps_cell[i][j] / self.lammps_tilt[k][k]) > 0.5
                and self.pbc[k]
                for i, j, k in FLIP_ORDER
            ]
        )
        for iteri, (i, j, k) in enumerate(FLIP_ORDER):
            if self.flip[iteri]:
                change = self.lammps_cell[k][k]
                change *= np.sign(self.lammps_cell[i][j])
                self.lammps_cell[i][j] -= change

    def get_lammps_prism(self):
        """Return into lammps coordination system rotated cell

        :returns: lammps cell
        :rtype: np.array

        """
        return self.lammps_cell[(0, 1, 2, 1, 2, 2), (0, 1, 2, 0, 0, 1)]

    # !TODO: detect flip in lammps
    def update_cell(self, lammps_cell):
        """Rotate new lammps cell into ase coordinate system

        :param lammps_cell: new lammps cell recieved after executing lammps
        :returns: ase cell
        :rtype: np.array

        """
        self.lammps_cell = lammps_cell
        self.lammps_tilt = self.lammps_cell.copy()

        # reverse flip
        for iteri, (i, j, k) in enumerate(FLIP_ORDER):
            if self.flip[iteri]:
                change = self.lammps_cell[k][k]
                change *= np.sign(self.lammps_cell[i][j])
                self.lammps_tilt[i][j] -= change

        # try to detect potential flips in lammps
        # (lammps minimizes the cell-vector lenghts)
        new_ase_cell = np.dot(self.lammps_tilt, self.rot_mat.T)
        # assuming the cell changes are mostly isotropic
        new_vol = np.linalg.det(new_ase_cell)
        old_vol = np.linalg.det(self.ase_cell)
        test_residual = self.ase_cell.copy()
        test_residual *= (new_vol / old_vol) ** (1.0 / 3.0)
        test_residual -= new_ase_cell
        if any(
                np.linalg.norm(test_residual, axis=1)
                > 0.5 * np.linalg.norm(self.ase_cell, axis=1)
        ):
            print(
                "WARNING: Significant simulation cell changes from LAMMPS "
                + "detected.\n"
                + " " * 9
                + "Backtransformation to ASE might fail!"
            )
        return new_ase_cell

    def vector_to_lammps(self, vec, wrap=False):
        """Rotate vector from ase coordinate system to lammps one

        :param vec: to be rotated ase-vector
        :returns: lammps-vector
        :rtype: np.array

        """
        # !TODO: right eps-limit
        # lammps might not like atoms outside the cell
        if wrap:
            return wrap_positions(
                np.dot(vec, self.rot_mat),
                self.lammps_cell,
                pbc=self.pbc,
                eps=1e-18,
            )

        return np.dot(vec, self.rot_mat)

    def vector_to_ase(self, vec, wrap=False):
        """Rotate vector from lammps coordinate system to ase one

        :param vec: to be rotated lammps-vector
        :param wrap: was vector wrapped into 'lammps_cell'
        :returns: ase-vector
        :rtype: np.array

        """
        if wrap:
            # trying to reverse wraping across periodic boundaries in the
            # lammps coordination-system:
            # translate: expresses lammps-coordinate system in the rotate, but
            #            without tilt removed system
            # fractional: vector in tilted system
            translate = np.linalg.solve(
                self.lammps_tilt.T, self.lammps_cell.T
            ).T
            fractional = np.linalg.solve(self.lammps_tilt.T, vec.T).T

            # !TODO: make somehow nicer
            # !TODO: handle extreme tilted cells
            for ifrac in fractional:
                for zyx in reversed(range(3)):
                    if ifrac[zyx] >= 1.0 and self.pbc[zyx]:
                        ifrac -= translate[zyx]
                    elif ifrac[zyx] < 0.0 and self.pbc[zyx]:
                        ifrac += translate[zyx]
            vec = np.dot(fractional, self.lammps_tilt)

        return np.dot(vec, self.rot_mat.T)

    def is_skewed(self):
        """Test if a lammps cell is not tetragonal

        :returns: bool
        :rtype: bool

        """
        cell_sq = self.lammps_cell ** 2
        return (
            np.sum(np.tril(cell_sq, -1)) / np.sum(np.diag(cell_sq))
            > self.tolerance
        )
