#!/usr/bin/env python3

'''
TODO: Description Needed
This modules tests aims calculator on different machines - Hawk, Isambard and archer.
The dimensions from 0 to 3 test the different aims setting as described in
carmm.run.get_aims_calculator

'''

def test_run_aims():
    from carmm.run.aims_path import set_aims_command

    for hpc in ['hawk', 'isambard', 'archer']:
        set_aims_command(hpc)

    from ase.calculators.aims import Aims
    from carmm.run.aims_calculator import get_aims_and_sockets_calculator

    for state in range(4):
        #fhi_calc = get_aims_calculator(state)
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(state, verbose=True)

        # Assertion test that the correct calculators are being set
        assert(type(sockets_calc.calc) == Aims)

test_run_aims()
