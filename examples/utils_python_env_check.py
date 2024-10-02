"""
Test script for python_env_check

This is a utility class for checking Python/ASE environments, and so wouldn't normally be used standalone.
"""

def test_ase_env_check():

    import ase
    from carmm.utils.python_env_check import ase_env_check

    assert ase_env_check(ase.__version__) == True

    # Checks if defaults are working
    if ase.__version__ >= '3.23.0':
        assert ase_env_check() == True
    else:
        assert ase_env_check() == False

def test_python_env_check():

    import sys
    from carmm.utils.python_env_check import python_env_check

    assert python_env_check(sys.version_info.minor) == True
    assert python_env_check(int(sys.version_info.minor)+1) == False

    # Checks if defaults are working
    if sys.version_info.minor >= 7:
        assert python_env_check() == True
    else:
        assert python_env_check() == False

def test_is_env_python():

    import sys
    from carmm.utils.python_env_check import is_env_python
    
    assert is_env_python(str(sys.version_info.minor)) == True
    assert is_env_python(str(int(sys.version_info.minor)+1)) == False

    assert is_env_python(int(sys.version_info.minor)) == True
    assert is_env_python(int(sys.version_info.minor)+1) == False

    test_value = (int(sys.version_info.major), int(sys.version_info.minor))
    assert is_env_python(test_value) == True
    test_value = (test_value[0]+1, test_value[1])
    assert is_env_python(test_value) == False

    try: 
        is_env_python(None) 
    except TypeError: 
        pass

test_ase_env_check()
test_python_env_check()
test_is_env_python()

