#!/usr/bin/env python3

'''
TODO: Description Needed

'''

def test_run_aims():
    from carmm.run.aims_path import set_aims_command

    set_aims_command('hawk')
    set_aims_command('isambard')
    set_aims_command('archer')

    from carmm.run.aims_calculator import get_aims_calculator

    calc = get_aims_calculator("gas")
    calc = get_aims_calculator("periodic")

    # TODO: Add an assertion test

test_run_aims()

