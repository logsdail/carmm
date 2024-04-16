import numpy as np
from carmm.build.neb.geodesic_utils import get_scaled_bond_dist_and_deriv, Morse, cost_func, cost_func_der
from ase.geometry import get_distances
class GeodesicInterpolator:
    def __init__(self, initial, final, path_len):

        from ase.constraints import FixAtoms

        self.initial = initial
        self.final = final
        self.images = [self.initial, self.final]
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

    def init_path(self):

        from ase.build.rotate import minimize_rotation_and_translation

        if self.path_len > 12:
            path_probe_len = self.path_len
        else:
            path_probe_len = 12

        for im_idx in range(path_probe_len):
            distances = [self.cart_rmsd(im1, im2) for im1, im2 in zip(self.images[1:], self.images)]
            mid_idx = np.argmax(distances)
            #print(f"{distances}, with {mid_idx} selected")

            new_mid = self.geodesic_midpoint(self.images[mid_idx], self.images[mid_idx + 1])
            self.images.insert(mid_idx + 1, new_mid)

            minimize_rotation_and_translation(self.images[mid_idx], self.images[mid_idx + 1])
            for im in range(len(self.images) - 1):
                minimize_rotation_and_translation(self.images[im], self.images[im + 1])

        # Sequentially remove images until correct number found
        while len(self.images) > self.path_len:
            dists = [self.cart_rmsd(im1, im2) for im1, im2 in zip(self.images[2:], self.images)]
            min_rmsd = np.argmin(dists)

            #print(f"Removing image {min_rmsd + 1}")
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

        from ase.neighborlist import natural_cutoffs, NeighborList

        cutoffs = np.array(natural_cutoffs(image)) * scale_cutoffs
        nl_init = NeighborList(cutoffs, self_interaction=False, bothways=False)
        nl_init.update(image)

        conn_mat = nl_init.get_connectivity_matrix()

        # filter for constraints
        if constraint:
            bond_idx = [tuple(sorted(bond)) for bond in conn_mat.keys() if
                        set() == set(self.fixed).intersection(set(bond))]
        else:
            bond_idx = [tuple(sorted(bond)) for bond in conn_mat.keys()]

        return bond_idx

    def geodesic_midpoint(self, image1, image2):

        from scipy.optimize import least_squares

        bonds_im1 = self.get_atoms_bonds(image1)
        bonds_im2 = self.get_atoms_bonds(image2)
        morse = Morse(alpha=0.7, re=1.5, beta=0.001)
        friction = 0.1 / np.sqrt(len(image1))

        comb_bonds = list(set(bonds_im1) | set(bonds_im2))
        new_bonds = []

        while True:

            d_min = np.inf

            comb_bonds = list(set(comb_bonds + new_bonds))

            init_wij, init_dwij = get_scaled_bond_dist_and_deriv(image1, comb_bonds, morse=morse)
            final_wij, final_dwij = get_scaled_bond_dist_and_deriv(image2, comb_bonds, morse=morse)

            av_wij = (init_wij + final_wij) / 2

            trial_atoms = image1.copy()

            for coef in [0.02, 0.98]:

                trial_atoms.positions = (image1.positions * coef) + ((1 - coef) * image2.positions)
                trial_atoms[self.movers].rattle(stdev=0.001)

                init_trial_pos0 = trial_atoms.positions.copy().ravel()
                init_trial_pos = trial_atoms.positions.copy().ravel()

                result = least_squares(cost_func,
                                       init_trial_pos,
                                       cost_func_der,
                                       args=(trial_atoms, comb_bonds, av_wij, friction, morse),
                                       ftol=0.02,
                                       gtol=0.02, verbose=0)

                opt_trial_pos = result.x.reshape((len(trial_atoms), 3))
                trial_atoms.positions = opt_trial_pos

                trial_bonds = self.get_atoms_bonds(trial_atoms)

                new_bonds = list(set(trial_bonds) - set(comb_bonds))

                if new_bonds:
                    #print(f"New bond detected: {new_bonds}. Trying again.")
                    break

                segment_len = self.compute_displacement_segment(image1, trial_atoms, image2, comb_bonds, morse=morse)

                left_len = self.cart_rmsd(image1, trial_atoms)
                right_len = self.cart_rmsd(image2, trial_atoms)
                width = max(left_len, right_len)

                dist = width + segment_len
                if dist < d_min:
                    d_min = dist
                    output_trial_atoms = trial_atoms.copy()

            else:
                break

        return output_trial_atoms

    def compute_displacement_segment(self, start_im, mid_im, end_im, bond_idx, friction=1e-3, morse=None):

        init_wij, init_dwij = get_scaled_bond_dist_and_deriv(start_im, bond_idx, morse=morse)
        mid_wij, mid_dwij = get_scaled_bond_dist_and_deriv(mid_im, bond_idx, morse=morse)
        final_wij, final_dwij = get_scaled_bond_dist_and_deriv(end_im, bond_idx, morse=morse)

        left_vec = mid_wij - init_wij
        right_vec = final_wij - mid_wij

        length_segment = np.sum(np.linalg.norm(left_vec, axis=1)) + np.sum(np.linalg.norm(right_vec, axis=1))

        return length_segment