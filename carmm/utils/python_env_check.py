"""
Created on Fri 12/07/24

@author Firas Assa

This function checks if the environment version entered as an argument is true.
"""

import sys
import warnings


def deprecation(message):
    # sub-function of python_env_check()
   warnings.warn(message, DeprecationWarning, stacklevel=2)


def python_env_check(version=None):
    while version is not None:
        if sys.version_info.minor >= version:
            return True
        else:
            return False
    else:
        print('\n' + sys.version)
        if sys.version_info.minor == 7:
           deprecation("Warning! Python Interpreter is too old, will be deprecated soon.")
        elif sys.version_info.minor > 7:
           print("Python Interpreter version is up-to-date.")
        else:
           assert sys.version_info >= (3, 6)


def is_env_python_minor(version):
    #sub-function of is_env_python()
    if sys.version_info.minor == version:
        print("1")
    else:
        print("0")


def is_env_python_major_and_minor(version):
    # sub-function of is_env_python()
    if sys.version_info.major == version[0] and sys.version_info.minor == version[1]:
        print("1")
    else:
        print("0")


def is_env_python(version):
    if isinstance(version, str):
        version = int(version)
        is_env_python_minor(version)
    elif isinstance(version, int):
        is_env_python_minor(version)
    elif isinstance(version, tuple):
        is_env_python_major_and_minor(version)
    else:
        raise TypeError('Invalid data type: enter only str, int or tuple')