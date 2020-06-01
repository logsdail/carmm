#!/usr/bin/env python3

'''
TODO: Description Needed
This modules tests aims calculator on different machines - Hawk, Isambard and archer.
The dimensions from 0 to 3 test the different aims setting on as described in
carmm.run.get_aims_calculator

'''

def test_run_aims():
    from carmm.run.aims_path import set_aims_command

    for hpc in ['hawk', 'isambard', 'archer']:
        set_aims_command(hpc)

    from carmm.run.aims_calculator import get_aims_calculator

    for state in range(4):
        calc = get_aims_calculator(state)

    # TODO: Add an actual assertion test

test_run_aims()
