from ase.geometry import get_distances, get_distances_derivatives
import numpy as np

def get_bond_dist_and_deriv(atoms, bond_list, input_pos=None, deriv_only=False, constraints=None):
    """
    Calculates distances between bonds in bond list and the derivatives of corresponding bond vectors
    with respect to cartesian coordinates of all atoms

    Parameters
    ----------
    atoms : ASE Atoms Object
        Atom information in ASE format
    bond_list : list
        List of tuples containing all considered bonds
    input_pos : ndarray
        1D array of positions, replacing positions of atoms (dim: len(atoms)*3)
    deriv_only : bool
        Returns only derivatives
    constraints : list
        List of fixed atomic positions

    Returns
    -------
    distances_arr : ndarray
        Array of distances between atoms specified in bond_list (dim: len(bond_list))
    derivatives_arr : ndarray
        Array of derivatives of bond vectors by cartesian coordinates of all atoms (dim: [len(bond_list), len(atoms), 3])
    """

    if input_pos is not None:
        save_positions = atoms.positions.copy()
        atoms.positions = input_pos.reshape(-1, 3)

    bond_array = np.array(bond_list)
    distances_arr = np.zeros(shape=(len(bond_list)))
    derivatives_arr = np.zeros(shape=(len(bond_list), len(atoms), 3))

    distances, vects = atoms.get_all_distances(mic=True), atoms.get_all_distances(mic=True, vector=True)
    distances_arr = np.array([ distances[idx1, idx2] for idx1, idx2 in zip(bond_array[:,0], bond_array[:,1]) ])

    bond_vects = np.array([ vects[idx1, idx2, :] for idx1, idx2 in zip(bond_array[:,0], bond_array[:,1]) ])
    derivatives = get_distances_derivatives(bond_vects, cell=atoms.get_cell(), pbc=atoms.pbc)

    for idx, bond in enumerate(bond_list):
        derivatives_arr[idx, bond[0], :] = derivatives[idx, 0, :]
        derivatives_arr[idx, bond[1], :] = derivatives[idx, 1, :]

    if constraints is not None:
        derivatives_arr[:, constraints, :] = np.array([0, 0, 0])
        derivatives_arr[:, constraints, :] = np.array([0, 0, 0])

    if input_pos is not None:
        atoms.positions = save_positions

    if deriv_only:
        return derivatives_arr

    return distances_arr, derivatives_arr


def get_scaled_bond_dist_and_deriv(atoms, bond_list, morse=None, input_pos=None, constraints=None):
    """
    Calculates distances between bonds in bond list and the derivatives of corresponding bond vectors
    with respect to cartesian coordinates of all atoms, scaled by an approximate Morse potential.

    Parameters
    ----------
    atoms : ASE Atoms Object
        Atom information in ASE format
    bond_list : list
        List of tuples containing all considered bonds
    morse : Class Morse
        Morse potential, which scales both distances and derivatives depending on pairwise distance and species
    input_pos : ndarray
        1D array of positions, replacing positions of atoms (dim: len(atoms)*3)
    constraints : list
        List of fixed atomic positions

    Returns
    -------
    distances_arr : ndarray
        Array of distances between atoms specified in bond_list, scaled by Morse potential (dim: len(bond_list))
    derivatives_arr : ndarray
        Array of derivatives of bond vectors by cartesian coordinates of all atoms, scaled by Morse potential
        (dim: [len(bond_list) * len(atoms), 3]).
    """

    distance_arr, derivative_arr = get_bond_dist_and_deriv(atoms, bond_list, input_pos=input_pos, constraints=constraints)

    if morse is None:
        # If not specified, just use default Morse parameters
        morse = Morse()

    morse.get_re(atoms, bond_list)

    dist_func = np.frompyfunc(morse.potential_en, 2, 1)
    deriv_func = np.frompyfunc(morse.potential_der, 2, 1)

    scaled_distance_arr = dist_func(distance_arr, morse.re_list, out=np.zeros(distance_arr.shape), casting='unsafe')
    scale_deriv = deriv_func(distance_arr, morse.re_list, out=np.zeros(distance_arr.shape), casting='unsafe')

    scaled_derivative_arr = np.einsum('ijk, i -> ijk', derivative_arr, scale_deriv)

    return scaled_distance_arr, scaled_derivative_arr.reshape(len(bond_list), -1)

class Morse:
    def __init__(self, alpha=0.7, beta=0.01):
        self.alpha = alpha
        self.beta = beta

    def potential_en(self, dist, re):
        ratio = dist / re
        val1 = np.exp(self.alpha * (1 - ratio))
        val2 = self.beta / ratio

        return val1 + val2

    def potential_der(self, dist, re):
        ratio = dist / re
        val1 = np.exp(self.alpha * (1 - ratio))
        val2 = self.beta / ratio
        deriv_scale = -(val1 * self.alpha / re) - (val2 / dist)

        return deriv_scale

    def get_re(self, atoms, bond_list):

        from ase.data import covalent_radii

        self.re_list = []
        for bond in bond_list:
            bond_idx1, bond_idx2 = bond[0], bond[1]
            atomic_number_1, atomic_number_2 = atoms.numbers[bond_idx1], atoms.numbers[bond_idx2]

            self.re_list.append(covalent_radii[atomic_number_1] + covalent_radii[atomic_number_2])

        self.re_list = np.array(self.re_list)

def cost_func(flat_positions, atoms, bond_list, av_wij, friction, morse, constraints):
    """
    Calculatos cost-function for minimisation of geodesic distance for newly generated image
    between two images

    Parameters
    ----------
    flat_positions : ndarray
        1D array of positions, replacing positions of atoms (dim: len(atoms)*3)
    atoms : ASE Atoms Object
        Atom information in ASE format
    bond_list : list
        List of tuples containing all considered bonds
    av_wij : ndarray
        Averaged internal coordinate distance between the preceding and succeeding image of the optimised
        image
    friction : float64
        Scaling factor
    morse : Class Morse
        Morse potential, which scales both distances and derivatives depending on pairwise distance and species
    constraints : list
        List of fixed atomic positions

    Returns
    -------
    output_wij:
        Cost function for least squares problem (dim: len(bond_list)+len(atoms))
    """

    wij, dwij = get_scaled_bond_dist_and_deriv(atoms, bond_list, input_pos=flat_positions, morse=morse, constraints=constraints)
    wij = (wij - av_wij)

    output_wij = np.concatenate([wij, (flat_positions - (atoms.get_positions().ravel())) * friction])

    return output_wij

def cost_func_der(flat_positions, atoms, bond_list, av_wij, friction, morse, constraints):
    """
    Calculatos gradient of cost-function for minimisation of geodesic distance for newly generated image
    between two images

    Parameters
    ----------
    flat_positions : ndarray
        1D array of positions, replacing positions of atoms (dim: len(atoms)*3)
    atoms : ASE Atoms Object
        Atom information in ASE format
    bond_list : list
        List of tuples containing all considered bonds
    av_wij : ndarray
        Averaged internal coordinate distance between the preceding and succeeding image of the optimised
        image
    friction : float64
        Scaling factor
    morse : Class Morse
        Morse potential, which scales both distances and derivatives depending on pairwise distance and species
    constraints : list
        List of fixed atomic positions

    Returns
    -------
    np.vstack([dwij, friction_scale]) : ndarray
        Cost function gradient for least squares problem (dim:  [len(bond_list) * len(atoms) + len(atoms), 3])
    """

    wij, dwij = get_scaled_bond_dist_and_deriv(atoms, bond_list, input_pos=flat_positions, morse=morse, constraints=constraints)

    friction_scale = np.identity(flat_positions.size) * friction

    return np.vstack([dwij, friction_scale])

