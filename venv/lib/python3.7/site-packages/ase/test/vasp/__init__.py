from ase.test import require
import os
import unittest


def installed():
    require('vasp')
    for env in ['VASP_COMMAND', 'VASP_SCRIPT']:
        if os.getenv(env):
            break
    else:
        raise unittest.SkipTest('Neither VASP_COMMAND nor VASP_SCRIPT defined')
    return True


def installed2():
    # Check if env variables exist for Vasp2
    require('vasp')

    for env in ['VASP_COMMAND', 'VASP_SCRIPT', 'ASE_VASP_COMMAND']:
        if os.getenv(env):
            break
    else:
        raise unittest.SkipTest('Neither ASE_VASP_COMMAND, VASP_COMMAND nor VASP_SCRIPT defined')

    return True
