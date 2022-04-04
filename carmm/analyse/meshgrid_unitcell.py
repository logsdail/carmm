class UnitCell:
    '''
    Object which stores and defines the unit cell parameters (number of points and dimensions) and the
    mesh grid of the underlying unit cell. Primarily used to reduce the number of variables which
    need to be passed to associated meshgrid functions. Currently only defined for orthogonal unit cells.
    '''
    def __init__(self):

        self.nx, self.ny, self.nz = 0, 0, 0
        self.X, self.Y, self.Z = 0, 0, 0
        self.xx, self.yy, self.zz = 0, 0, 0

        self.dim = 0

        self.box_means = 0
        self.box_mean_indices = 0

    def define_unit_cell(self, nx, ny, nz, x, y, z):
        '''
        Defines a unit cell of (nx, ny, nz) points of (x, y, z) length.
        Args:
            n_x, n_y, n_z: int
                Number of points along each cartesian axis.
            x, y, z: float
                Length of each cartesian axis.
        '''

        import numpy as np

        # Define the number of points and the dimensions of each cartesian axis.
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dim = np.array([x, y, z])

        # Define the linespace of points spanned by each cartesian axis.
        self.X = np.linspace(0, self.dim[0], self.nx)
        self.Y = np.linspace(0, self.dim[1], self.ny)
        self.Z = np.linspace(0, self.dim[2], self.nz)

        # Define meshgrid of the unit cell.
        self.xx, self.yy, self.zz = np.meshgrid(self.X, self.Y, self.Z, indexing='xy')