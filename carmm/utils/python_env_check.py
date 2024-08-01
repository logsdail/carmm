"""
Created on Fri 12/07/24

@author Firas Assa
"""

def deprecation(message, stacklevel=2):
    """
    A general purpose deprecation warning for use when raising future error warnings

    Args:
        message: str
            Text to display with warning
        stacklevel: int
            Sets the number of lines to be unwound when raising the warning.

    Returns:
        Nothing
            

    """

    import warnings

    # sub-function of python_env_check()
    warnings.warn(message, DeprecationWarning, stacklevel=stacklevel)

def ase_env_check(version=None):
    """
    A general purpose check for the ASE version being used, to filter code routes

    Args:
        version: str
            Text representation of the version number for ASE

    Returns:
        Boolean: true if version is compatible, false if not

    """

    import ase

    if version is None:
        version = '3.23.0' # Base requirement for functionality

    if ase.__version__ >= version:
        return True
    else:
        deprecation("Warning! Your current ASE version is: " + ase.__version__ + ", and depreciated. v3.23 or higher is recommended.")
        return False

def python_env_check(minorversion=None):
    """
    A general purpose check for the Python minor version being used, to filter code routes

    Args:
        minorversion: int
            Minor version number of Python to check. Currently assumes major version is 3.

    Returns:
        Boolean: true if version is compatible, false if not

    """
    import sys

    if minorversion is None:
        minorversion = 7

    if sys.version_info.minor >= minorversion:
        return True
    else:
        deprecation("Warning! Your current Python Interpreter is: " + sys.version + ", and will be deprecated soon.")
        return False

def _is_env_python_minor(version):
    """
    A sub-function of is_env_python, checking minor version exact match

    Args:
        version: int
            Minor version number of Python to check. Currently assumes major version is 3.

    Returns:
        Boolean: true if version is identical, false if not

    """

    import sys

    if sys.version_info.minor == version:
        return True
    else:
        return False

def _is_env_python_major_and_minor(version):
    """
    A sub-function of is_env_python, checking major and minor version exact match

    Args:
        version: tuple, int
            Major and minor version number of Python to check.

    Returns:
        Boolean: true if version is identical, false if not

    """


    import sys

    if sys.version_info.major == version[0] and sys.version_info.minor == version[1]:
        return True
    else:
        return False

def is_env_python(version):
    """
    A helper function to check exact match of Python version

    Args:
        version: int, str, tuple
            Version number of Python to check. Can handle various formats to aid compatibility

    Returns:
        Boolean: true if version is identical, false if not

    """

    if isinstance(version, str):
        version = int(version)
        return _is_env_python_minor(version)
    elif isinstance(version, int):
        return _is_env_python_minor(version)
    elif isinstance(version, tuple):
        return _is_env_python_major_and_minor(version)
    else:
        raise TypeError('Invalid data type: enter only str, int or tuple')
