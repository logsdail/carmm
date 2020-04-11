import sys
import numpy as np
from pymysql import connect
from pymysql.err import ProgrammingError
from copy import deepcopy

from ase.db.sqlite import SQLite3Database
from ase.db.sqlite import init_statements
from ase.db.sqlite import VERSION
from ase.db.postgresql import remove_nan_and_inf, insert_nan_and_inf
import ase.io.jsonio
import json


class Connection(object):
    """
    Wrapper for the MySQL connection

    Arguments
    =========
    host: str
        Hostname. For a local database this is localhost.
    user: str
        Username.
    passwd: str
        Password
    db_name: str
        Name of the database
    port: int
        Port
    binary_prefix: bool
        MySQL checks if an argument can be interpreted as a UTF-8 string. This
        check fails for binary values. Binary values need to have _binary
        prefix in MySQL. By setting this to True, the prefix is automatically
        added for binary values.
    """

    def __init__(self, host=None, user=None, passwd=None, port=3306,
                 db_name=None, binary_prefix=False):
        self.con = connect(host=host, user=user, passwd=passwd, db=db_name,
                           binary_prefix=binary_prefix)

    def cursor(self):
        return MySQLCursor(self.con.cursor())

    def commit(self):
        self.con.commit()

    def close(self):
        self.con.close()


class MySQLCursor(object):
    """
    Wrapper for the MySQL cursor. The most important task performed by this
    class is to translate SQLite queries to MySQL. Translation is needed
    because ASE DB uses some field names that are reserved words in MySQL.
    Thus, these has to mapped onto other field names.
    """
    sql_replace = [
        (' key TEXT', ' attribute_key TEXT'),
        ('(key TEXT', '(attribute_key TEXT'),
        ('SELECT key FROM', 'SELECT attribute_key FROM'),
        ('?', '%s'),
        (' keys ', ' attribute_keys '),
        (' key=', ' attribute_key='),
        ('table.key', 'table.attribute_key'),
        (' IF NOT EXISTS', '')
    ]

    def __init__(self, cur):
        self.cur = cur

    def execute(self, sql, params=None):

        # Replace external table key -> attribute_key
        for substibution in self.sql_replace:
            sql = sql.replace(substibution[0], substibution[1])

        if params is None:
            params = ()

        self.cur.execute(sql, params)

    def fetchone(self):
        return self.cur.fetchone()

    def fetchall(self):
        return self.cur.fetchall()

    def _replace_nan_inf_kvp(self, values):
        for item in values:
            if not np.isfinite(item[1]):
                item[1] = sys.float_info.max / 2
        return values

    def executemany(self, sql, values):
        if 'number_key_values' in sql:
            values = self._replace_nan_inf_kvp(values)

        for substibution in self.sql_replace:
            sql = sql.replace(substibution[0], substibution[1])
        self.cur.executemany(sql, values)


class MySQLDatabase(SQLite3Database):
    """
    ASE interface to a MySQL database (via pymysql package).

    Arguments
    ==========
    url: str
        URL to the database. It should have the form
        mysql://username:password@host:port/database_name.
        Example URL with the following credentials
            username: john
            password: johnspasswd
            host: localhost (i.e. server is running locally)
            database: johns_calculations
            port: 3306
        mysql://john:johnspasswd@localhost:3306/johns_calculations
    create_indices: bool
        Carried over from parent class. Currently indices are not
        created for MySQL, as TEXT fields cannot be hashed by MySQL.
    use_lock_file: bool
        See SQLite
    serial: bool
        See SQLite
    """
    type = 'mysql'
    default = 'DEFAULT'

    def __init__(self, url=None, create_indices=True,
                 use_lock_file=False, serial=False):
        super(MySQLDatabase, self).__init__(
            url, create_indices, use_lock_file, serial)

        self.host = None
        self.username = None
        self.passwd = None
        self.db_name = None
        self.port = 3306
        self._parse_url(url)

    def _parse_url(self, url):
        """
        Parse the URL
        """
        url = url.replace('mysql://', '')
        url = url.replace('mariadb://', '')

        splitted = url.split(':', 1)
        self.username = splitted[0]

        splitted = splitted[1].split('@')
        self.passwd = splitted[0]

        splitted = splitted[1].split('/')
        host_and_port = splitted[0].split(':')
        self.host = host_and_port[0]
        self.port = int(host_and_port[1])
        self.db_name = splitted[1]

    def _connect(self):
        return Connection(host=self.host, user=self.username,
                          passwd=self.passwd, db_name=self.db_name,
                          port=self.port, binary_prefix=True)

    def _initialize(self, con):
        if self.initialized:
            return

        cur = con.cursor()

        information_exists = True
        self._metadata = {}
        try:
            cur.execute("SELECT 1 FROM information")
        except ProgrammingError:
            information_exists = False

        if not information_exists:
            # We need to initialize the DB
            # MySQL require that id is explicitly set as primary key
            # in the systems table
            init_statements_cpy = deepcopy(init_statements)
            init_statements_cpy[0] = init_statements_cpy[0][:-1] + \
                ', PRIMARY KEY(id))'

            statements = schema_update(init_statements_cpy)
            for statement in statements:
                cur.execute(statement)
            con.commit()
            self.version = VERSION
        else:
            cur.execute('select * from information')

            for name, value in cur.fetchall():
                if name == 'version':
                    self.version = int(value)
                elif name == 'metadata':
                    self._metadata = json.loads(value)

        self.initialized = True

    def blob(self, array):
        if array is None:
            return None
        return super(MySQLDatabase, self).blob(array).tobytes()

    def get_offset_string(self, offset, limit=None):
        sql = ''
        if not limit:
            # mysql does not allow for setting limit to -1 so
            # instead we set a large number
            sql += '\nLIMIT 10000000000'
        sql += '\nOFFSET {0}'.format(offset)
        return sql

    def get_last_id(self, cur):
        cur.execute('select max(id) as ID from systems')
        last_id = cur.fetchone()[0]
        return last_id

    def create_select_statement(self, keys, cmps,
                                sort=None, order=None, sort_table=None,
                                what='systems.*'):
        sql, value = super(MySQLDatabase, self).create_select_statement(
            keys, cmps, sort, order, sort_table, what)

        for subst in MySQLCursor.sql_replace:
            sql = sql.replace(subst[0], subst[1])
        return sql, value

    def encode(self, obj, binary=False):
        return ase.io.jsonio.encode(remove_nan_and_inf(obj))

    def decode(self, obj, lazy=False):
        return insert_nan_and_inf(ase.io.jsonio.decode(obj))


def schema_update(statements):
    for i, statement in enumerate(statements):
        for a, b in [('REAL', 'DOUBLE'),
                     ('INTEGER PRIMARY KEY AUTOINCREMENT',
                      'INT NOT NULL AUTO_INCREMENT')]:
            statements[i] = statement.replace(a, b)

    # MySQL does not support UNIQUE constraint on TEXT
    # need to use VARCHAR. The unique_id is generated with
    # randint(16**31, 16**32-1) so it will contain 32
    # hex-characters
    statements[0] = statements[0].replace('TEXT UNIQUE', 'VARCHAR(32) UNIQUE')

    # keys is a reserved word in MySQL redefine this table name to
    # attribute_keys
    statements[2] = statements[2].replace('keys', 'attribute_keys')

    txt2jsonb = ['calculator_parameters', 'key_value_pairs']

    for column in txt2jsonb:
        statements[0] = statements[0].replace(
            '{} TEXT,'.format(column),
            '{} JSON,'.format(column))

    statements[0] = statements[0].replace('data BLOB,', 'data JSON,')

    tab_with_key_field = ['attribute_keys', 'number_key_values',
                          'text_key_values']

    # key is a reserved word in MySQL redefine this to attribute_key
    for i, statement in enumerate(statements):
        for tab in tab_with_key_field:
            if tab in statement:
                statements[i] = statement.replace(
                    'key TEXT', 'attribute_key TEXT')
    return statements
