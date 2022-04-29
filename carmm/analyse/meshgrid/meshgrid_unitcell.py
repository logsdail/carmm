class MeshObject:
    """
    Object which stores and defines the unit cell parameters (number of points and dimensions) and the
    mesh grid of the underlying unit cell. Primarily used to reduce the number of variables which
    need to be passed to associated meshgrid functions.
    """

    def __init__(self, cell_dims, nx=50, ny=50, nz=50, pbc=[0,0,0]):

        import numpy as np
        from ase.geometry import Cell, cellpar_to_cell

        if np.shape(cell_dims)==(6,):
            cellpar = cellpar_to_cell(cell_dims)
            self.Cell = Cell(cellpar)
        elif np.shape(cell_dims)==(3,3,3):
            self.Cell = Cell(cell_dims)
        else:
            raise Exception("Specify cell parameters as a (6,) (x,y,z,alpha,beta,gamma) or (3x3x3) array.")

        self.nx, self.ny, self.nz = nx, ny, nz
        self.pbc = pbc

        self.X = np.linspace(0, 1, self.nx)
        self.Y = np.linspace(0, 1, self.ny)
        self.Z = np.linspace(0, 1, self.nz)

        # Define meshgrid in fractional coordinates.
        self.frac_xx, self.frac_yy, self.frac_zz = np.meshgrid(self.X, self.Y, self.Z, indexing='xy')

        # Project from fractional coordinates to cartesian.
        self.xx = np.dot(self.frac_xx, self.Cell.array)
        self.yy = np.dot(self.frac_yy, self.Cell.array)
        self.zz = np.dot(self.frac_zz, self.Cell.array)