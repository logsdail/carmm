"""
Test if we can find vasp2 using get_calculator()
"""

from ase.test.vasp import installed2 as installed
from ase.calculators.calculator import get_calculator_class

assert installed()

get_calculator_class('vasp2')
