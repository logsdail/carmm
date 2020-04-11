from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError


class LinearCombinationCalculator(Calculator):
    """LinearCombinationCalculator for weighted summation of multiple calculators.
    """

    def __init__(self, calcs, weights, atoms=None):
        """Implementation of sum of calculators.

        calcs: list
            List of an arbitrary number of :mod:`ase.calculators` objects.
        weights: list of float
            Weights for each calculator in the list.
        atoms: Atoms object
            Optional :class:`~ase.Atoms` object to which the calculator will be attached.
        """

        super().__init__(atoms=atoms)

        if len(calcs) == 0:
            raise ValueError('The value of the calcs must be a list of Calculators')

        for calc in calcs:
            if not isinstance(calc, Calculator):
                raise ValueError('All the calculators should be inherited form the ase\'s Calculator class')

        common_properties = set.intersection(*(set(calc.implemented_properties) for calc in calcs))
        self.implemented_properties = list(common_properties)

        if not self.implemented_properties:
            raise PropertyNotImplementedError('There are no common property implemented for the potentials!')

        if len(weights) != len(calcs):
            raise ValueError('The length of the weights must be the same as the number of calculators!')

        self.calcs = calcs
        self.weights = weights

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super().calculate(atoms, properties, system_changes)

        if not set(properties).issubset(self.implemented_properties):
            raise PropertyNotImplementedError('Some of the requested property is not in the '
                                              'list of supported properties ({})'.format(self.implemented_properties))

        for w, calc in zip(self.weights, self.calcs):
            calc.calculate(atoms, properties, system_changes)

            for k in properties:
                if k not in self.results:
                    self.results[k] = w * calc.results[k]
                else:
                    self.results[k] += w * calc.results[k]

    def reset(self):
        """Clear all previous results recursively from all fo the calculators."""
        super().reset()

        for calc in self.calcs:
            calc.reset()

    def __str__(self):
        calculators = ', '.join(calc.__class__.__name__ for calc in self.calcs)
        return '{}({})'.format(self.__class__.__name__, calculators)


class MixedCalculator(LinearCombinationCalculator):
    """
    Mixing of two calculators with different weights

    H = weight1 * H1 + weight2 * H2

    Has functionality to get the energy contributions from each calculator

    Parameters
    ----------
    calc1 : ASE-calculator
    calc2 : ASE-calculator
    weight1 : float
        weight for calculator 1
    weight2 : float
        weight for calculator 2
    """

    def __init__(self, calc1, calc2, weight1, weight2):
        super().__init__([calc1, calc2], [weight1, weight2])

    def set_weights(self, w1, w2):
        self.weights[0] = w1
        self.weights[1] = w2

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super().calculate(atoms, properties, system_changes)
        if 'energy' in properties:
            energy1 = self.calcs[0].get_property('energy', atoms)
            energy2 = self.calcs[1].get_property('energy', atoms)
            self.results['energy_contributions'] = (energy1, energy2)

    def get_energy_contributions(self, atoms=None):
        """ Return the potential energy from calc1 and calc2 respectively """
        self.calculate(properties=['energy'], atoms=atoms)
        return self.results['energy_contributions']


class SumCalculator(LinearCombinationCalculator):
    """SumCalculator for combining multiple calculators.

    This calculator can be used when there are different calculators for the different chemical environment or
    for example during delta leaning. It works with a list of arbitrary calculators and evaluates them in sequence
    when it is required.
    The supported properties are the intersection of the implemented properties in each calculator.
    """

    def __init__(self, calcs, atoms=None):
        """Implementation of sum of calculators.

        calcs: list
            List of an arbitrary number of :mod:`ase.calculators` objects.
        atoms: Atoms object
            Optional :class:`~ase.Atoms` object to which the calculator will be attached.
        """

        weights = [1.] * len(calcs)
        super().__init__(calcs, weights, atoms)


class AverageCalculator(LinearCombinationCalculator):
    """AverageCalculator for equal summation of multiple calculators (for thermodynamic purposes)..
    """

    def __init__(self, calcs, atoms=None):
        """Implementation of average of calculators.

        calcs: list
            List of an arbitrary number of :mod:`ase.calculators` objects.
        atoms: Atoms object
            Optional :class:`~ase.Atoms` object to which the calculator will be attached.
        """

        n = len(calcs)

        if n == 0:
            raise ValueError('The value of the calcs must be a list of Calculators')

        weights = [1 / n] * n
        super().__init__(calcs, weights, atoms)
