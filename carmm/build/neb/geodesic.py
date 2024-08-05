import numpy as np
from carmm.build.neb.geodesic_utils import get_scaled_bond_dist_and_deriv, Morse, cost_func, cost_func_der
from ase.geometry import get_distances

class GeodesicInterpolator:
    """
    Method for producing NEB pathways using the geodesic interpolation method of
    Zhu et al. [1]. Adapted to ASE objects, and operates with periodic boundary
    conditions and FixedAtom constraints of ASE.

    Example
    -------
    from ase.io import read
    from carmm.build.neb.geodesic import GeodesicInterpolator
    from ase.constraints import FixAtoms
    from ase.mep import NEB

    # Example below produces a pathway with 6 interpolating images
    initial = read("initial.xyz")
    final = read("final.xyz")

    # If constraints are desired, specify as normal
    constraints = Fixatoms([0, 1])
    initial.set_constraints(constraints)
    final.set_constraints(constraints)

    # Interpolate images by iteratively generating images between the midpoint
    # images along the pathway
    geodesic = GeodesicInterpolator(initial, final, path_len=8)
    geodesic.init_path()

    # Smooth the pathway by performing a sweeping algorithm which reduces the
    # overall geodesic length of the pathway
    geodesic.sweep_iterative(sweeperiter=20)

    # Define the NEB pathway as normal
    neb = NEB(geodesic.images)

    Citations
    ---------
    [1] Xiaolei Zhu et al.; Geodesic interpolation for reaction pathways. J. Chem. Phys.
        28 April 2019; 150 (16): 164103.
    """
    def __init__(self, initial, final, path_len):

        from ase.constraints import FixAtoms
        from ase.build.rotate import minimize_rotation_and_translation

        self.initial = initial
        self.final = final
        self.images = [self.initial, self.final]
        minimize_rotation_and_translation(self.images[0], self.images[1])

        self.cell = self.initial.get_cell()
        self.pbc = self.initial.pbc
        self.path_len = path_len

        # Check constraints - only honour fixed atoms
        constraints = self.initial._get_constraints()

        self.fixed = []
        for constraint in constraints:
            if type(constraint) == FixAtoms:
                self.fixed = [a for a in constraint.index]

        self.movers = list(
            set(np.arange(len(self.initial))) - set(self.fixed) | set(self.fixed) - set(np.arange(len(self.initial))))

        self.bond_list = []
        self.update_bond_list()

        self.initialized_internals = False

    def init_path(self):
        """
        Initialises the pathway by iteratively adding images at the midpoint between two
        existing images, such that the geodesic length between the resulting three
        points is minimized.

        Algorithm structured as follows:
        1) For 18 given iterations, find the maximum RMSD distance between two images
        2) Select two said images and find the mid-point image
        3) Align images using a Kabsch algorithm to minimise rotations and translations
        """

        from ase.build.rotate import minimize_rotation_and_translation
        from ase.visualize import view

        if self.path_len > 20:
            path_probe_len = self.path_len
        else:
            path_probe_len = 20

        for im_idx in range(path_probe_len-2):
            distances = [self.cart_rmsd(im2, im1) for im1, im2 in zip(self.images[1:], self.images)]
            mid_idx = np.argmax(distances)
            print(f"Inserting new image with distance, {distances}, with {mid_idx} selected")

            new_mid = self.geodesic_midpoint(self.images[mid_idx], self.images[mid_idx + 1])
            #new_mid = self.geodesic_midpoint(self.images[mid_idx], self.images[mid_idx + 1], image_selection=im_idx%2)
            self.update_bond_list()
            self.images.insert(mid_idx + 1, new_mid)

            for im in np.arange(1, len(self.images) - 1):
                minimize_rotation_and_translation(self.images[im - 1], self.images[im])
            self.images = list(reversed(self.images))
            for im in np.arange(1, len(self.images) - 1):
                minimize_rotation_and_translation(self.images[im - 1], self.images[im])
            self.images = list(reversed(self.images))

        # Sequentially remove images until correct number found
        while len(self.images) > self.path_len:
            dists = [self.cart_rmsd(im1, im2) for im1, im2 in zip(self.images[2:], self.images)]
            min_rmsd = np.argmin(dists)

            print(f"Removing image {min_rmsd + 1}")
            del self.images[min_rmsd + 1]

        # Ensure constraints are honoured
        for im_idx, image in enumerate(self.images):
            self.images[im_idx].positions[self.fixed] = self.initial.positions[self.fixed]

    def cart_rmsd(self, atoms1, atoms2):

        dist_cart = \
        get_distances(atoms1.positions[self.movers], atoms2.positions[self.movers], pbc=self.pbc, cell=self.cell)[1]
        rmsd = np.mean(np.diag(dist_cart))

        return rmsd

    def get_atoms_bonds(self, image, scale_cutoffs=1., constraint=True):

        from ase.neighborlist import natural_cutoffs, NeighborList, NewPrimitiveNeighborList

        cutoffs = np.array(natural_cutoffs(image)) * scale_cutoffs
        nl_init = NeighborList(cutoffs, self_interaction=False, bothways=False, primitive=NewPrimitiveNeighborList)
        nl_init.update(image)

        conn_mat = nl_init.get_connectivity_matrix()

        # filter for constraints
        if constraint:
            '''
            If constraints are present, only include bonds either between non-constrained
            atoms or a constrained atom and non-constrained atom.
            
            The cost function gradient for constrained atoms is set to 0, in principle, 
            only allowing non-constrained atom to vary position in least squares 
            geodesic interpolation.
            '''
            bond_idx = [tuple(sorted(bond)) for bond in conn_mat.keys() if
                        2 > len(list(set(self.fixed).intersection(set(bond))))
                        and bond[0] != bond[1]]
        else:
            bond_idx = [tuple(sorted(bond)) for bond in conn_mat.keys()]

        return bond_idx

    def smooth_iteration(self, im_idx1, im_idx2, xref=None):

        from ase.build.rotate import minimize_rotation_and_translation
        from scipy.optimize import least_squares
        from ase.visualize import view

        morse = Morse(alpha=1.7, beta=0.1)
        friction = 0.001

        #self.update_internal_coordinates(morse)

        flat_pos_0 = self.images[im_idx1].positions.ravel().copy()

        result = least_squares(self.cost_func_2atoms,
                               flat_pos_0,
                               self.cost_func_der_2atoms,
                               args=(im_idx1, im_idx2, friction, morse, xref),
                               ftol=0.02,
                               gtol=0.02, verbose=0,
                               loss='soft_l1')
        #print(f"Optimality: {self.segment_length, self.optimality}")
        #print(flat_pos_0.reshape((len(self.images[im_idx1]),3)) - result.x.reshape((len(self.images[im_idx1]), 3)))
        self.images[im_idx1].positions = result.x.reshape((len(self.images[im_idx1]), 3))

        #self.update_bond_list()
        #minimize_rotation_and_translation(self.images[im_idx2], self.images[im_idx1])

    def update_bond_list(self):

        for idx in range(len(self.images)):
            bonds = self.get_atoms_bonds(self.images[idx], scale_cutoffs=3.0)
            #bonds = get_bond_list(self.images[idx].positions.ravel(), threshold=3)[0]
            for bond in bonds:
                self.bond_list.append(bond)

        self.bond_list = list(set(self.bond_list))
        self.initialized_internals = False

    def update_internal_coordinates(self, morse, update_idx = None):

        if not self.initialized_internals:
            print("Clearing Wij data")
            self.wij_list = np.zeros((len(self.images), len(self.bond_list)))
            self.dwij_list = np.zeros( (len(self.images), len(self.bond_list), len(self.images[0]) * 3) )
            self.mid_wij_list = np.zeros((len(self.images) - 1, len(self.bond_list)))
            self.mid_dwij_list = np.zeros( (len(self.images) - 1, len(self.bond_list), len(self.images[0]) * 3) )
            self.mid_images = self.images[1:]
            self.initialized_internals = True

        if update_idx is None:
            update_images = self.images
            update_idx = np.arange(len(self.images))
        else:
            update_images = [self.images[image] for image in update_idx]

        for idx, image in zip(update_idx, update_images):
            wij, dwij = get_scaled_bond_dist_and_deriv(image, self.bond_list, morse=morse, constraints=self.fixed)
            self.wij_list[idx, :] = wij
            self.dwij_list[idx, :] = dwij

        for idx, image in zip(update_idx, update_images):
            if idx == len(self.images) - 1:
                continue
            self.mid_images[idx] = image.copy()
            self.mid_images[idx].positions = (self.images[idx].positions + self.images[idx+1].positions)/2

        for idx, image in zip(update_idx, update_images):
            if idx == len(self.images) - 1:
                continue
            wij, dwij = get_scaled_bond_dist_and_deriv(self.mid_images[idx], self.bond_list, morse=morse, constraints=self.fixed)
            self.mid_wij_list[idx, :] = wij
            self.mid_dwij_list[idx, :] = dwij

    def sweep_iterative(self, sweeperiter=20):
        """
        Adjusts the iterpolating images to minimize the geodesic length over a
        given segment in the pathway. Sweeps are alternatively performed from
        reactant to product and product to reactant

        Parameters
        ----------
        sweeperiter : int
            Number of sweeping iterations performed
        -------
        """

        from ase.build.rotate import minimize_rotation_and_translation
        from ase.visualize import view

        ''' Ensures that the path is correctly oriented on exit '''
        if sweeperiter%2 == 1:
            sweeperiter = sweeperiter + 1

        for iter in range(sweeperiter):

            # Forward sweep
            for idx in np.arange(1, len(self.images)-1):
                print(f"Forward sweep iter: {iter}, idx: {idx}")

                xmid = (self.images[idx-1].positions + self.images[idx+1].positions)/2
                self.smooth_iteration(idx, idx+1, xref=xmid)

                # Backward sweep
                self.images = list(reversed(self.images))

        for im in np.arange(1, len(self.images) - 1):
            minimize_rotation_and_translation(self.images[im-1], self.images[im])
        self.images = list(reversed(self.images))
        for im in np.arange(1, len(self.images) - 1):
            minimize_rotation_and_translation(self.images[im-1], self.images[im])
        self.images = list(reversed(self.images))

        #for im in np.arange(1, len(self.images)-1 ):
        #    minimize_rotation_and_translation(self.images[im + 1],self.images[im])

    def geodesic_midpoint(self, image1, image2, image_selection=None):

        from ase.build.rotate import minimize_rotation_and_translation
        from scipy.optimize import least_squares
        from ase.visualize import view

        morse = Morse(alpha=0.7, beta=0.01)
        friction = 0.1 / np.sqrt(len(image1))

        output_trial_atoms = None
        new_bonds = []
        bonds_im1 = self.get_atoms_bonds(image1, scale_cutoffs=3.0)
        bonds_im2 = self.get_atoms_bonds(image2, scale_cutoffs=3.0)
        #bonds_im1 = get_bond_list(image1.positions.ravel(), threshold=4 + 1, enforce=set(new_bonds))[0]
        #bonds_im2 = get_bond_list(image1.positions.ravel(), threshold=4 + 1, enforce=set(new_bonds))[0]
        comb_bonds = list(set(bonds_im1) | set(bonds_im2))

        while output_trial_atoms is None:

            d_min = np.inf
            left_dist = np.inf
            right_dist = np.inf
            new_bonds_detected = False

            comb_bonds = list(set(comb_bonds + new_bonds))

            init_wij, init_dwij = get_scaled_bond_dist_and_deriv(image1, comb_bonds, morse=morse)
            final_wij, final_dwij = get_scaled_bond_dist_and_deriv(image2, comb_bonds, morse=morse)
            av_wij = (init_wij + final_wij) / 2

            trial_atoms = image1.copy()

            for coef in [0.02, 0.98]:
            #for random in range(10):
                trial_atoms.positions = (image1.positions * coef) + ((1 - coef) * image2.positions)
                trial_atoms.rattle(stdev=0.01)
                #print(random)
                #print(rand)
                #if random%2 == 0:
                #    input_pos = image1.positions + (0.05 * np.random.random_sample((len(image2),3)))
                #else:
                #    input_pos = image2.positions + (0.05 * np.random.random_sample((len(image2),3)))

                #trial_atoms[self.movers].positions += 0.1 * np.random.random_sample((len(self.movers),3))
                #trial_atoms.positions = input_pos

                init_trial_pos = trial_atoms.positions.ravel().copy()

                result = least_squares(cost_func,
                                       init_trial_pos,
                                       cost_func_der,
                                       args=(trial_atoms, comb_bonds, av_wij, friction, morse, self.fixed),
                                       ftol=0.02,
                                       gtol=0.02, verbose=0)
                                       #loss="soft_l1")

                opt_trial_pos = result.x.reshape((len(trial_atoms), 3))
                trial_atoms.positions = opt_trial_pos

                trial_bonds = self.get_atoms_bonds(trial_atoms, scale_cutoffs=3.0)
                #trial_bonds = get_bond_list(trial_atoms.positions.ravel(), threshold=4 + 1, enforce=set(new_bonds))[0]

                new_bonds = list(set(trial_bonds) - set(comb_bonds))

                if new_bonds:
                    print(f"New bond detected: {new_bonds}. Trying again.")
                    new_bonds_detected = True
                    break

                old_images = self.images
                self.images = [image1, trial_atoms, image2]
                self.sweep_iterative(sweeperiter=1)
                trial_atoms = self.images[1]
                self.images = old_images

                segment_len = self.compute_displacement_segment(image1, trial_atoms, image2, comb_bonds, morse=morse)
                left_len = self.cart_rmsd(image1, trial_atoms)
                right_len = self.cart_rmsd(image2, trial_atoms)
                #centrality = abs((abs(left_len - right_len))/(left_len + right_len) - 0.5)
                width = max(left_len, right_len)
                dist = width + segment_len

                if dist < d_min:
                    print(f"Coef {coef} accepted with dist {dist}")
                    d_min = dist
                    output_trial_atoms = trial_atoms.copy()

        return output_trial_atoms

    def compute_displacement_all(self, start_im, mid_im, end_im, bond_idx, friction=1e-3, morse=None):

        mid_1_pos = (start_im.positions.ravel() + mid_im.positions.ravel())/2
        mid_2_pos = (end_im.positions.ravel() + mid_im.positions.ravel())/2

        init_wij, init_dwij = get_scaled_bond_dist_and_deriv(start_im, bond_idx, morse=morse)
        init_mid_wij, init_mid_dwij = get_scaled_bond_dist_and_deriv(start_im, bond_idx, input_pos=mid_1_pos, morse=morse)
        mid_wij, mid_dwij = get_scaled_bond_dist_and_deriv(mid_im, bond_idx, morse=morse)
        final_wij, final_dwij = get_scaled_bond_dist_and_deriv(end_im, bond_idx, morse=morse)
        final_mid_wij, final_mid_dwij = get_scaled_bond_dist_and_deriv(end_im, bond_idx, input_pos=mid_1_pos, morse=morse)

        seg_1 = np.linalg.norm(init_wij - init_mid_wij) + np.linalg.norm(init_mid_wij - mid_wij)
        seg_2 = np.linalg.norm(mid_wij - final_mid_wij) + np.linalg.norm(final_mid_wij - final_wij)
        #seg_3 = np.linalg.norm(mid_wij - init_mid_wij) + np.linalg.norm(init_mid_wij - mid_wij)

        #left_vec = (mid_wij - init_wij)
        #right_vec = (final_wij - mid_wij)

        #length_segment = np.sum(np.linalg.norm(left_vec, axis=1)) + np.sum(np.linalg.norm(right_vec, axis=1))
        #length_segment = np.sum(np.linalg.norm(left_vec) + np.linalg.norm(right_vec))
        length_segment = seg_1 + seg_2

        return length_segment

    def compute_displacement_segment(self, start_im, mid_im, end_im, bond_idx, friction=1e-3, morse=None):

        mid_1_pos = (start_im.positions.ravel() + mid_im.positions.ravel())/2
        mid_2_pos = (end_im.positions.ravel() + mid_im.positions.ravel())/2

        init_wij, init_dwij = get_scaled_bond_dist_and_deriv(start_im, bond_idx, morse=morse)
        init_mid_wij, init_mid_dwij = get_scaled_bond_dist_and_deriv(start_im, bond_idx, input_pos=mid_1_pos, morse=morse)
        mid_wij, mid_dwij = get_scaled_bond_dist_and_deriv(mid_im, bond_idx, morse=morse)
        final_wij, final_dwij = get_scaled_bond_dist_and_deriv(end_im, bond_idx, morse=morse)
        final_mid_wij, final_mid_dwij = get_scaled_bond_dist_and_deriv(end_im, bond_idx, input_pos=mid_1_pos, morse=morse)

        seg_1 = np.linalg.norm(init_wij - init_mid_wij) + np.linalg.norm(init_mid_wij - mid_wij)
        seg_2 = np.linalg.norm(mid_wij - final_mid_wij) + np.linalg.norm(final_mid_wij - final_wij)
        #seg_3 = np.linalg.norm(mid_wij - init_mid_wij) + np.linalg.norm(init_mid_wij - mid_wij)

        #left_vec = (mid_wij - init_wij)
        #right_vec = (final_wij - mid_wij)

        #length_segment = np.sum(np.linalg.norm(left_vec, axis=1)) + np.sum(np.linalg.norm(right_vec, axis=1))
        #length_segment = np.sum(np.linalg.norm(left_vec) + np.linalg.norm(right_vec))
        length_segment = seg_1 + seg_2

        return length_segment

    def cost_func_2atoms(self, flat_positions, im_idx1, im_idx2, friction, morse, xmid):

        self.images[im_idx1].positions = np.reshape(flat_positions, (len(self.images[im_idx1]), 3))
        self.update_internal_coordinates(morse, [im_idx1-1, im_idx1, im_idx2])

        #wij_1 = np.sqrt((self.mid_wij_list[im_idx1] - self.wij_list[im_idx1])**2)
        #wij_2 = np.sqrt((self.wij_list[im_idx2] - self.mid_wij_list[im_idx1])**2)
        wij_1 = (self.mid_wij_list[im_idx1 - 1] - self.wij_list[im_idx1 - 1])
        wij_2 = (self.mid_wij_list[im_idx1] - self.wij_list[im_idx1])
        wij_3 = (self.wij_list[im_idx1] - self.mid_wij_list[im_idx1 - 1])
        wij_4 = (self.wij_list[im_idx2] - self.mid_wij_list[im_idx1])

        self.segment_length = np.sum(np.linalg.norm(wij_1) + np.linalg.norm(wij_2))

        #output_wij = wij_1 + wij_2
        output_wij = np.concatenate([wij_1, wij_2, wij_3, wij_4, ((flat_positions - xmid.ravel()) * friction)])

        self.two_atom_wij = output_wij.copy()

        return self.two_atom_wij

    def cost_func_der_2atoms(self, flat_positions, im_idx1, im_idx2, friction, morse, xmid):

        grad = np.zeros((4 * len(self.bond_list) + 3 * (im_idx2 - im_idx1) * len(self.images[im_idx1]), 3 * len(self.images[im_idx1])))

        mid1 = self.mid_dwij_list[im_idx1 - 1] / 2
        mid2 = self.mid_dwij_list[im_idx1] / 2
        grad[(1 * len(self.bond_list)):(2 * len(self.bond_list)), :] = mid2 - self.dwij_list[im_idx1]
        grad[:(1 * len(self.bond_list)), :] = mid1
        grad[(3 * len(self.bond_list)):(4 * len(self.bond_list)), :] = -mid2
        grad[(2 * len(self.bond_list)):(3 * len(self.bond_list)), :] = self.dwij_list[im_idx1] - mid1

        for idx in range((im_idx2 - im_idx1) * 3 * len(self.images[im_idx1])):
            grad[4 * len(self.bond_list) + idx, idx] = friction

        self.two_atom_dwij = grad
        #self.optimality = np.linalg.norm(np.einsum('i,i...', self.two_atom_wij, self.two_atom_dwij), ord=np.inf)

        return self.two_atom_dwij