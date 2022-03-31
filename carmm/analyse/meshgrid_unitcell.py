class Unit_Cell:
    '''
    Object which stores and defines the unit cell parameters (number of points and dimensions) and the
    mesh grid of the underlying unit cell. Primarily used to reduce the number of variables which
    need to be passed to associated meshgrid functions. Currently only defined for orthogonal unit cells.
    '''

    def __init__(self):

        self.nx, self.ny, self.nz = 0, 0, 0
        self.X, self.Y, self.Z = 0, 0, 0
        self.xx, self.yy, self.zz = 0, 0, 0

        self.ind_list_x, self.ind_list_y, self.ind_list_z = 0, 0, 0
        self.n_part_box = 0
        self.part_x, self.part_y, self.part_z = 0, 0, 0
        self.part_nx, self.part_ny, self.part_nz = 0, 0, 0
        self.part_dx, self.part_dy, self.part_dz = 0, 0, 0

        self.dim = 0

        self.box_means = 0
        self.box_mean_indices = 0
        self.parallel = False

    def define_unit_cell(self, x, y, z, nx=50, ny=50, nz=50, parallel=False, part_x=0, part_y=0, part_z=0):
        '''
        Defines a unit cell of (nx, ny, nz) points of (x, y, z) length. Also defines the partition box indices
        used for splitting up the simulation cell into smaller units.
        Args:
            n_x, n_y, n_z: int
                Number of points along each cartesian axis.
            x, y, z: float
                Length of each cartesian axis.
        '''

        import numpy as np

        # Define the number of points and the dimensions of each cartesian axis.
        if (nx < 10 or ny < 10 or nz < 10):
            raise ValueError('Number of points along one axis too low. Please increase above 9.')
        if (nx < 1 or ny < 1 or nz < 1):
            raise ValueError('Number of points along one or more axis negative.')

        self.nx, self.ny, self.nz = nx, ny, nz
        self.dim = np.array([x, y, z])

        # Define the linespace of points spanned by each cartesian axis.
        self.X = np.linspace(0, self.dim[0], self.nx)
        self.Y = np.linspace(0, self.dim[1], self.ny)
        self.Z = np.linspace(0, self.dim[2], self.nz)

        # Define meshgrid of the unit cell.
        self.xx, self.yy, self.zz = np.meshgrid(self.X, self.Y, self.Z, indexing='xy')

        self.parallel = parallel

        if parallel:
            print(f"Parallelisation specified. Looking for mpi4py...")
            print(f"mpi4py detected. Defining partition boxes...")

        if (part_x == 0 or part_y == 0 or part_z == 0):
            print(f"Automatically searching for number of partition boxes...")
            self.find_partition_nboxes()
            self.define_partition_boxes(self.part_x, self.part_y, self.part_z)
        else:
            self.define_partition_boxes(part_x, part_y, part_z)

    def find_partition_nboxes(self):

        import numpy as np

        x_factors = np.array([int(x) for x in np.arange(np.round((self.nx + 1) / 2)) if self.nx % x == 0])
        y_factors = np.array([int(y) for y in np.arange(np.round((self.ny + 1) / 2)) if self.ny % y == 0])
        z_factors = np.array([int(z) for z in np.arange(np.round((self.nz + 1) / 2)) if self.nz % z == 0])

        # Try and get the number of points per box as close to 10 as possible...
        if self.nx < 20:
            self.part_x = int(self.nx / x_factors[np.argmax(x_factors)])
        else:
            self.part_x = int(self.nx / x_factors[np.argmin(np.absolute(x_factors - 10))])

        if self.ny < 20:
            self.part_y = int(self.ny / y_factors[np.argmax(y_factors)])
        else:
            self.part_y = int(self.ny / y_factors[np.argmin(np.absolute(y_factors - 10))])

        if self.nz < 20:
            self.part_z = int(self.nz / z_factors[np.argmax(z_factors)])
        else:
            self.part_z = int(self.nz / z_factors[np.argmin(np.absolute(z_factors - 10))])

        if (self.part_x == self.nx or self.part_y == self.ny or self.part_z == self.nz):
            print("Warning: One or more of your axis are defined only using one partition box.\n")
            print("Code efficiency can be improved setting nx/ny/nz which can be divided into even groups of points.")

    def define_partition_boxes(self, n_box_x, n_box_y, n_box_z):

        # Define the number of partition boxes along each axis.
        self.part_x, self.part_y, self.part_z = n_box_x, n_box_y, n_box_z

        # Define the length of each partition box axis.
        self.part_dx, self.part_dy, self.part_dz = self.nx / self.dim[0], self.ny / self.dim[1], self.nz / self.dim[2]

        # Define the total number of partition boxes.
        self.n_part_box = self.part_x * self.part_y * self.part_z

        # Gather start and end indices of points used in each partition box.
        ind_list_x = []
        ind_list_y = []
        ind_list_z = []
        self.ind_list_x = self.partition_indices(self.part_x, self.nx, ind_list_x)
        self.ind_list_y = self.partition_indices(self.part_y, self.ny, ind_list_y)
        self.ind_list_z = self.partition_indices(self.part_z, self.nz, ind_list_z)

        # Define the number of points along each axis of the partition boxes.
        self.part_nx = self.ind_list_x[1][0] - self.ind_list_x[0][0]
        self.part_ny = self.ind_list_y[1][0] - self.ind_list_y[0][0]
        self.part_nz = self.ind_list_z[1][0] - self.ind_list_z[0][0]

        # Calculate the mean positions of all partition boxes.
        self.calculate_part_box_means(self)

    def calculate_part_box_means(self):

        import numpy as np

        self.box_means = []
        self.box_mean_indices = []
        for x in range(self.part_x):
            for y in range(self.part_y):
                for z in range(self.part_z):
                    x_mean = np.mean(self.xx[self.ind_list_x[x][0]:self.ind_list_x[x][1],
                                     self.ind_list_y[y][0]:self.ind_list_y[y][1],
                                     self.ind_list_z[z][0]:self.ind_list_z[z][1]])
                    y_mean = np.mean(self.yy[self.ind_list_x[x][0]:self.ind_list_x[x][1],
                                     self.ind_list_y[y][0]:self.ind_list_y[y][1],
                                     self.ind_list_z[z][0]:self.ind_list_z[z][1]])
                    z_mean = np.mean(self.zz[self.ind_list_x[x][0]:self.ind_list_x[x][1],
                                     self.ind_list_y[y][0]:self.ind_list_y[y][1],
                                     self.ind_list_z[z][0]:self.ind_list_z[z][1]])
                    self.box_mean_indices.append([x, y, z])
                    self.box_means.append([x_mean, y_mean, z_mean])

    def partition_indices(self, partition_axis, npoints, ind_list):

        # Define the start and end axes of the partition boxes.
        p_per_box = npoints / partition_axis

        for npart in range(partition_axis):
            min_partind_ax = int(npart * p_per_box)
            max_partind_ax = int(((npart + 1) * p_per_box))

            ind_list.append([min_partind_ax, max_partind_ax])

        return ind_list
