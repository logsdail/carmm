from ase.test import must_raise
import unittest
try:
    import pymysql
    _ = pymysql  # Supress unused import warning
except ImportError:
    raise unittest.SkipTest('No MySQL module')


from ase.db import connect
from ase import Atoms
from ase.calculators.emt import EMT
from ase.build import molecule
import os

ON_CI_SERVER = 'CI_PROJECT_DIR' in os.environ.keys()

if ON_CI_SERVER:
    URL = 'mysql://root:ase@mysql:3306/testase_mysql'
    # HOST = 'mysql'
    # USER = 'root'
    # PASSWD = 'ase'
    # DB_NAME = 'testase_mysql'
else:
    URL = os.environ.get('MYSQL_DB_URL')
    # HOST = os.environ.get('MYSQL_HOST', None)
    # USER = os.environ.get('MYSQL_USER', None)
    # PASSWD = os.environ.get('MYSQL_PASSWD', None)
    # DB_NAME = os.environ.get('MYSQL_DB_NAME', None)

if URL is None:
    raise unittest.SkipTest('Not on GitLab CI server. To run this test '
                            'host, username, password and database name '
                            'must be in the environment variables '
                            'MYSQL_HOST, MYSQL_USER, MYSQL_PASSWD and '
                            'MYSQL_DB_NAME, respectively.')


def full_db_name():
    return URL


def test_connect():
    db = connect(full_db_name())
    db.delete([row.id for row in db.select()])


def test_write_read():
    db = connect(full_db_name())

    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])
    uid = db.write(co, tag=1, type='molecule')

    co_db = db.get(id=uid)
    atoms_db = co_db.toatoms()

    assert len(atoms_db) == 2
    assert atoms_db[0].symbol == co[0].symbol
    assert atoms_db[1].symbol == co[1].symbol
    assert co_db.tag == 1
    assert co_db.type == 'molecule'


def test_write_read_with_calculator():
    db = connect(full_db_name())

    h2o = molecule('H2O')
    calc = EMT(dummy_param=2.4)
    h2o.set_calculator(calc)

    uid = db.write(h2o)

    h2o_db = db.get(id=uid).toatoms(attach_calculator=True)

    calc_db = h2o_db.get_calculator()
    assert calc_db.parameters['dummy_param'] == 2.4

    # Check that get_atoms function works
    db.get_atoms(H=2)


def test_update():
    db = connect(full_db_name())

    h2o = molecule('H2O')

    uid = db.write(h2o, type='molecule')
    db.update(id=uid, type='oxide')

    atoms_type = db.get(id=uid).type

    assert atoms_type == 'oxide'


def test_delete():
    db = connect(full_db_name())

    h2o = molecule('H2O')
    uid = db.write(h2o, type='molecule')

    # Make sure that we can get the value
    db.get(id=uid)
    db.delete([uid])

    with must_raise(KeyError):
        db.get(id=uid)


def test_read_write_bool_key_value_pair():
    db = connect(full_db_name())
    h2o = molecule('H2O')

    # Make sure we can read and write boolean key value pairs
    uid = db.write(h2o, is_water=True, is_solid=False)
    row = db.get(id=uid)
    assert row.is_water
    assert not row.is_solid

test_connect()
test_write_read()
test_write_read_with_calculator()
test_update()
test_delete()
test_read_write_bool_key_value_pair()
