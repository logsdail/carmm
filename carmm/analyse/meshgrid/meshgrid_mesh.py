class Mesh:
    """
    Object which stores and defines the unit cell parameters (number of points and dimensions) and the
    mesh grid of the underlying unit cell. Primarily used to reduce the number of variables which
    need to be passed to associated meshgrid functions.
    """

    def __init__(self, cell_dims, nx=50, ny=50, nz=50, pbc=[True,True,True], pbc_strict_mode=True):

        import numpy as np
        from ase.geometry import Cell, cellpar_to_cell


        if np.shape(cell_dims)==(6,):
            cellpar = cellpar_to_cell(cell_dims)
            self.Cell = Cell(cellpar)
        elif np.shape(cell_dims)==(3,3):
            self.Cell = Cell(cell_dims)
        else:
            raise Exception("Specify cell parameters as a (6,) (x,y,z,alpha,beta,gamma) or (3x3) array.")

        self.nx, self.ny, self.nz = nx, ny, nz

        self.pbc = np.zeros(3, bool)
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

        # Define the Cell parameters as (a,b,c,alpha,beta,gamma) array without calling
        # Cell.cellpar() method everytime.
        self.cellpar = self.Cell.cellpar()

        # Enforces Mesh and Atoms pbc match - if off, ignores mismatches.
        self.strict_mode = pbc_strict_mode

    def read_cube_file(self, cube_name):

        # Reads cube file and imports into a meshgrid object
        #
        # If cube file is structured Nx, Ny, Nz, and meshgrid xx, yy, zz.
        # To place data onto grid,
        #
        import numpy as np
        from ase.io.cube import read_cube

        cube_file = open(cube_name)
        cube_inp = read_cube(cube_file)

        c_shape = np.shape(cube_inp['data'])
        self.nx, self.ny, self.nz = c_shape[0], c_shape[1], c_shape[2]

        self.data = cube_inp.data
        atoms = cube_inp.atoms

        return atoms