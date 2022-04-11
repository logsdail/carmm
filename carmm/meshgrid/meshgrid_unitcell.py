class UnitCell:
    """
    Object which stores and defines the unit cell parameters (number of points and dimensions) and the
    mesh grid of the underlying unit cell. Primarily used to reduce the number of variables which
    need to be passed to associated meshgrid functions. Currently only defined for orthogonal unit cells.
    """

    def __init__(self):

        self.nx, self.ny, self.nz = 0, 0, 0
        self.X, self.Y, self.Z = 0, 0, 0
        self.xx, self.yy, self.zz = 0, 0, 0

        self.dim = 0

    def define_unit_cell(self, x, y, z, nx=50, ny=50, nz=50):
        """
        Defines a unit cell of (nx, ny, nz) points of (x, y, z) length. Also defines the partition box indices
        used for splitting up the simulation cell into smaller units.
        Args:
            n_x, n_y, n_z: int
                Number of points along each cartesian axis.
            x, y, z: float
                Length of each cartesian axis.
        """

        import numpy as np

        # Define the number of points and the dimensions of each cartesian axis.
        if nx < 10 or ny < 10 or nz < 10:
            raise ValueError('Number of points along one axis too low. Please increase above 9.')
        if nx < 1 or ny < 1 or nz < 1:
            raise ValueError('Number of points along one or more axis negative.')

        self.nx, self.ny, self.nz = nx, ny, nz
        self.dim = np.array([x, y, z])

        # Define the linespace of points spanned by each cartesian axis.
        self.X = np.linspace(0, self.dim[0], self.nx)
        self.Y = np.linspace(0, self.dim[1], self.ny)
        self.Z = np.linspace(0, self.dim[2], self.nz)

        # Define meshgrid of the unit cell.
        self.xx, self.yy, self.zz = np.meshgrid(self.X, self.Y, self.Z, indexing='xy')
