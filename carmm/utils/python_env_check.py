"""
Created on Fri 12/07/24

@author Firas Assa

deprecation() is a sub function of python_env_check() that outputs a DeprecationWarning

python_env_check() function has an optional argument version, which when is entered outputs True boolean value if python
interpreter version is >= current version. When no argument is entered, it checks whether python version > 3.7.

is_env_python() function checks if the environment version entered as an argument is true, by printing 1 or 0
"""

def deprecation(message):
    import warnings

    # sub-function of python_env_check()
    warnings.warn(message, DeprecationWarning, stacklevel=2)

def ase_env_check(version=None):
    import ase

    if version is None:
        version = '3.23.0' # Base requirement for functionality

    if ase.__version__ >= version:
        return True
    else:
        deprecation("Warning! Your current ASE version is: " + ase.__version__ + ", and depreciated. v3.23 or higher is recommended.")
        return False

def python_env_check(minorversion=None):
    import sys

    if minorversion is None:
        minorversion = 7

    if sys.version_info.minor >= minorversion:
        return True
    else:
        deprecation("Warning! Your current Python Interpreter is: " + sys.version + ", and will be deprecated soon.")
        return False

def is_env_python_minor(version):
    import sys

    #sub-function of is_env_python()
    if sys.version_info.minor == version:
        return True
    else:
        return False

def is_env_python_major_and_minor(version):
    import sys

    # sub-function of is_env_python()
    if sys.version_info.major == version[0] and sys.version_info.minor == version[1]:
        return True
    else:
        return False

def is_env_python(version):
    if isinstance(version, str):
        version = int(version)
        return is_env_python_minor(version)
    elif isinstance(version, int):
        return is_env_python_minor(version)
    elif isinstance(version, tuple):
        return is_env_python_major_and_minor(version)
    else:
        raise TypeError('Invalid data type: enter only str, int or tuple')
