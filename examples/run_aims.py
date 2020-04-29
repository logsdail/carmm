#!/usr/bin/env python3

'''
TODO: Description Needed

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

