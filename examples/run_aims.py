#!/usr/bin/env python3

'''
This modules tests aims calculator on different machines - Hawk, Isambard and archer.
The dimensions from 0 to 3 test the different aims setting as described in
carmm.run.get_aims_calculator
TODO: rework the assertion to actually test the aims_path output - include expected string values
'''

def test_run_aims():
    from carmm.run.aims_path import set_aims_command

    for hpc in ['hawk', 'isambard', 'archer2', 'young']:
        set_aims_command(hpc)

    import ase # Necessary to check the code version, as socket functionality has changed
    ase_major_version = int(ase.__version__.split(".")[0])
    ase_minor_version = int(ase.__version__.split(".")[1])

    from ase.calculators.aims import Aims
    from carmm.run.aims_calculator import get_aims_and_sockets_calculator

    for state in range(4):
        #fhi_calc = get_aims_calculator(state)
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(state, verbose=True)

        # Assertion test that the correct calculators are being set
        # ASE version 3.21 or earlier
        if ase_major_version <= 3 and ase_minor_version <= 21:
            assert (type(sockets_calc.calc) == Aims)
        else:
        # ASE Version 3.22 or later
            assert(type(sockets_calc.launch_client.calc) == Aims)

test_run_aims()
