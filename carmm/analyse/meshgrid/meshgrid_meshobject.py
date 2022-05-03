class Mesh:
    """
    Object which stores and defines the unit cell parameters (number of points and dimensions) and the
    mesh grid of the underlying unit cell. Primarily used to reduce the number of variables which
    need to be passed to associated meshgrid functions.
    """

    def __init__(self, cell_dims, nx=50, ny=50, nz=50, pbc=[0,0,0]):

        import numpy as np
        from ase.geometry import Cell, cellpar_to_cell, cell_to_cellpar

        if np.shape(cell_dims)==(6,):
            cellpar = cellpar_to_cell(cell_dims)
            self.Cell = Cell(cellpar)
        elif np.shape(cell_dims)==(3,3):
            self.Cell = Cell(cell_dims)
        else:
            raise Exception("Specify cell parameters as a (6,) (x,y,z,alpha,beta,gamma) or (3x3) array.")

        self.nx, self.ny, self.nz = nx, ny, nz
        self.pbc = pbc

        self.X = np.linspace(0, 1, self.nx)
        self.Y = np.linspace(0, 1, self.ny)
        self.Z = np.linspace(0, 1, self.nz)

        # Define meshgrid in fractional coordinates.
        self.frac_xx, self.frac_yy, self.frac_zz = np.meshgrid(self.X, self.Y, self.Z, indexing='xy')

        # Project from fractional coordinates to cartesian.
        stack_mesh = np.stack((self.frac_xx,self.frac_yy,self.frac_zz),axis=-1)
        stack_mesh = np.einsum('ji,abcj->iabc', self.Cell.array, stack_mesh)

        self.xx, self.yy, self.zz = stack_mesh[0], stack_mesh[1], stack_mesh[2]

        # Allows easy conversion from cart -> frac coordinates.
        self.inverse_cell_array = np.linalg.inv(self.Cell.array)

        self.cellpar = self.Cell.cellpar()
        self.x_max = np.dot(np.array([1,0,0]), self.Cell.array)[0] 
        self.y_max = np.dot(np.array([0,1,0]), self.Cell.array)[1]
        self.z_max = np.dot(np.array([0,0,1]), self.Cell.array)[2]
