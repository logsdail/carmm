from ase.geometry import get_distances, get_distances_derivatives
import numpy as np
def get_bond_dist_and_deriv(atoms, bond_list, input_pos=None, deriv_only=False):

    if input_pos is not None:
        atoms.positions = input_pos.reshape(-1, 3)

    distances_arr = np.zeros(shape=(len(atoms), len(atoms)))
    derivatives_arr = np.zeros(shape=(len(atoms), len(atoms), 3))

    for bond in bond_list:
        idx1, idx2 = bond[0], bond[1]

        dist = atoms.get_distance(idx1, idx2, mic=True)
        dist_vect = atoms.get_distance(idx1, idx2, mic=True, vector=True)
        deriv = get_distances_derivatives([dist_vect], cell=atoms.cell, pbc=atoms.pbc)[0][0][0]

        distances_arr[idx1, idx2] = dist
        derivatives_arr[idx1, idx2] = deriv

    if deriv_only:
        return derivatives_arr

    return distances_arr, derivatives_arr


def get_scaled_bond_dist_and_deriv(atoms, bond_list, morse=None, input_pos=None):
    distance_arr, derivative_arr = get_bond_dist_and_deriv(atoms, bond_list, input_pos=input_pos)

    if morse is None:
        # If not specified, just use default Morse parameters
        morse = Morse()

    dist_func = np.frompyfunc(morse.potential_en, 1, 1)
    deriv_func = np.frompyfunc(morse.potential_der, 1, 1)

    scaled_distance_arr = dist_func(distance_arr, where=distance_arr > 0, out=np.zeros(distance_arr.shape),
                                    casting='unsafe')

    scale_deriv = deriv_func(distance_arr, where=distance_arr > 0, out=np.zeros(distance_arr.shape), casting='unsafe')
    scaled_derivative_arr = np.einsum('ijk, ij -> ijk', derivative_arr, scale_deriv)

    return scaled_distance_arr, scaled_derivative_arr


class Morse:
    def __init__(self, re=1.5, alpha=0.7, beta=0.01):
        self.re = re
        self.alpha = alpha
        self.beta = beta

    def potential_en(self, dist):
        ratio = dist / self.re
        val1 = np.exp(self.alpha * (1 - ratio))
        val2 = self.beta / ratio

        return val1 + val2

    def potential_der(self, dist):
        ratio = dist / self.re
        val1 = np.exp(self.alpha * (1 - ratio))
        val2 = self.beta / ratio
        deriv_scale = -self.alpha / self.re * val1 - (val2 / dist)

        return deriv_scale


def cost_func(flat_positions, atoms, bond_list, av_wij, friction, morse):
    wij, dwij = get_scaled_bond_dist_and_deriv(atoms, bond_list, input_pos=flat_positions, morse=morse)
    wij0, dwij0 = get_scaled_bond_dist_and_deriv(atoms, bond_list, morse=morse)
    wij = wij - av_wij

    # Return as relevant bond list - otherwise least-squares problem blows up to natomsxnatomsxnatomsx3
    # for the Jacobian
    output_wij = np.zeros(len(bond_list))

    for idx, (i, j) in enumerate(bond_list):
        output_wij[idx] = wij[i, j]

    output_wij = np.concatenate([output_wij, (flat_positions - atoms.get_positions().ravel()) * friction])

    return output_wij


def cost_func_der(flat_positions, atoms, bond_list, av_wij, friction, morse):
    wij, dwij = get_scaled_bond_dist_and_deriv(atoms, bond_list, input_pos=flat_positions, morse=morse)

    # Return as relevant bond list - otherwise least-squares problem blows up to natomsxnatomsxnatomsx3
    # for the Jacobian
    output_dwij = np.zeros((len(bond_list), len(atoms), 3))

    for idx, (i, j) in enumerate(bond_list):
        output_dwij[idx, i] = dwij[i, j]
        output_dwij[idx, j] = -dwij[i, j]

    friction_scale = np.identity(flat_positions.size) * friction

    return np.vstack([output_dwij.reshape((len(bond_list), 3 * len(atoms))), friction_scale])
