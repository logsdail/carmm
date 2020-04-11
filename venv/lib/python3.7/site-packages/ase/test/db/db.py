import os
import time
from ase.test import cli
from ase.db import connect

cmd = """
ase -T build H | ase -T run emt -o testase.json &&
ase -T build H2O | ase -T run emt -o testase.json &&
ase -T build O2 | ase -T run emt -o testase.json &&
ase -T build H2 | ase -T run emt -f 0.02 -o testase.json &&
ase -T build O2 | ase -T run emt -f 0.02 -o testase.json &&
ase -T build -x fcc Cu | ase -T run emt -E 5,1 -o testase.json &&
ase -T db -v testase.json natoms=1,Cu=1 --delete --yes &&
ase -T db -v testase.json "H>0" -k hydro=1,abc=42,foo=bar &&
ase -T db -v testase.json "H>0" --delete-keys foo"""


def count(n, *args, **kwargs):
    m = len(list(con.select(columns=['id'], *args, **kwargs)))
    assert m == n, (m, n)


t0 = time.time()
for name in ['testase.json', 'testase.db', 'postgresql', 'mysql', 'mariadb']:
    if name == 'postgresql':
        if os.environ.get('POSTGRES_DB'):  # gitlab-ci
            name = 'postgresql://ase:ase@postgres:5432/testase'
        else:
            name = os.environ.get('ASE_TEST_POSTGRES_URL')
            if name is None:
                continue
    elif name == 'mysql':
        if os.environ.get('CI_PROJECT_DIR'):  # gitlab-ci
            name = 'mysql://root:ase@mysql:3306/testase_mysql'
        else:
            name = os.environ.get('MYSQL_DB_URL')

        if name is None:
            continue
    elif name == 'mariadb':
        if os.environ.get('CI_PROJECT_DIR'):  # gitlab-ci
            name = 'mariadb://root:ase@mariadb:3306/testase_mysql'
        else:
            name = os.environ.get('MYSQL_DB_URL')

        if name is None:
            continue

    con = connect(name)
    t1 = time.time()
    if 'postgres' in name or 'mysql' in name or 'mariadb' in name:
        con.delete([row.id for row in con.select()])

    cli(cmd.replace('testase.json', name))
    assert con.get_atoms(H=1)[0].magmom == 1
    count(5)
    count(3, 'hydro')
    count(0, 'foo')
    count(3, abc=42)
    count(3, 'abc')
    count(0, 'abc,foo')
    count(3, 'abc,hydro')
    count(0, foo='bar')
    count(1, formula='H2')
    count(1, formula='H2O')
    count(3, 'fmax<0.1')
    count(1, '0.5<mass<1.5')
    count(5, 'energy')

    id = con.reserve(abc=7)
    assert con[id].abc == 7

    for key in ['calculator', 'energy', 'abc', 'name', 'fmax']:
        count(6, sort=key)
        count(6, sort='-' + key)

    cli('ase -T gui --terminal {}@3'.format(name))

    con.delete([id])

    if name != 'testase.json':  # transfer between db formats
        if name == 'testase.db':
            from_db = 'testase.json'
            factor = 2
        else:
            from_db = 'testase.db'
            factor = 3

        cli('ase db {} --insert-into {}'.format(from_db, name))

        count(5 * factor)
        count(3 * factor, 'hydro')
        count(0, 'foo')
        count(3 * factor, abc=42)
        count(3 * factor, 'abc')
        count(0, 'abc,foo')
        count(3 * factor, 'abc,hydro')
        count(0, foo='bar')
        count(factor, formula='H2')
        count(factor, formula='H2O')
        count(3 * factor, 'fmax<0.1')
        count(factor, '0.5<mass<1.5')
        count(5 * factor, 'energy')

    t2 = time.time()

    print('----------------------------------')
    print('Finished test for {}'.format(name))
    print('runtime = {} sec'.format(t2 - t1))
    print('----------------------------------')


print('Total runtime = {} sec'.format(t2 - t0))
