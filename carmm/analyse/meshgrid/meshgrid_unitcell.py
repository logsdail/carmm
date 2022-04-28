class MeshObject:
    """
    Object which stores and defines the unit cell parameters (number of points and dimensions) and the
    mesh grid of the underlying unit cell. Primarily used to reduce the number of variables which
    need to be passed to associated meshgrid functions. Currently only defined for orthogonal unit cells.
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

        xdim, ydim, zdim = self.Cell.array[0][0], self.Cell.array[1][1], self.Cell.array[2][2]

        # Define the linespace of points spanned by each cartesian axis.
        self.X = np.linspace(0, xdim, self.nx)
        self.Y = np.linspace(0, ydim, self.ny)
        self.Z = np.linspace(0, zdim, self.nz)

        # Define meshgrid of the unit cell.
        self.xx, self.yy, self.zz = np.meshgrid(self.X, self.Y, self.Z, indexing='xy')

