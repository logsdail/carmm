import numpy as np
from ase import units
from ase.calculators.calculator import Calculator, all_changes


class SpringCalculator(Calculator):
    """
    Spring calculator corresponding to independent oscillators with a fixed
    spring constant.


    Energy for an atom is given as

    E = k / 2 * (r - r_0)**2

    where k is the spring constant and, r_0 the ideal positions.


    Parameters
    ----------
    ideal_positions : array
        array of the ideal crystal positions
    k : float
        spring constant in eV/Angstrom
    """
    implemented_properties = ['forces', 'energy']

    def __init__(self, ideal_positions, k):
        Calculator.__init__(self)
        self.ideal_positions = ideal_positions.copy()
        self.k = k

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        energy, forces = self.compute_energy_and_forces(atoms)
        self.results['energy'], self.results['forces'] = energy, forces

    def compute_energy_and_forces(self, atoms):
        disps = atoms.positions - self.ideal_positions
        forces = - self.k * disps
        energy = sum(self.k / 2.0 * np.linalg.norm(disps, axis=1)**2)
        return energy, forces

    def get_free_energy(self, T, method='classical'):
        """Get analytic vibrational free energy for the spring system.

        Parameters
        ----------
        T : float
            temperature (K)
        method : str
            method for free energy computation; 'classical' or 'QM'.
        """
        F = 0.0
        masses, counts = np.unique(self.atoms.get_masses(), return_counts=True)
        for m, c in zip(masses, counts):
            F += c * SpringCalculator.compute_Einstein_solid_free_energy(self.k, m, T, method)
        return F

    @staticmethod
    def compute_Einstein_solid_free_energy(k, m, T, method='classical'):
        """ Get free energy (per atom) for an Einstein crystal.

        Free energy of a Einstein solid given by classical (1) or QM (2)
        1.    F_E = 3NkbT log( hw/kbT )
        2.    F_E = 3NkbT log( 1-exp(hw/kbT) ) + zeropoint

        Parameters
        -----------
        k : float
            spring constant (eV/A^2)
        m : float
            mass (grams/mole or AMU)
        T : float
            temperature (K)
        method : str
            method for free energy computation, classical or QM.

        Returns
        --------
        float
            free energy of the Einstein crystal (eV/atom)
        """
        assert method in ['classical', 'QM']

        hbar = units._hbar * units.J  # eV/s
        m = m / units.kg              # mass kg
        k = k * units.m**2 / units.J  # spring constant J/m2
        omega = np.sqrt(k / m)        # angular frequency 1/s

        if method == 'classical':
            F_einstein = 3 * units.kB * T * np.log(hbar * omega / (units.kB * T))
        elif method == 'QM':
            log_factor = np.log(1.0 - np.exp(-hbar * omega / (units.kB * T)))
            F_einstein = 3 * units.kB * T * log_factor + 1.5 * hbar * omega

        return F_einstein
