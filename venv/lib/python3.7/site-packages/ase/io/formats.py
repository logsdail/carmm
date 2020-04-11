"""File formats.

This module implements the read(), iread() and write() functions in ase.io.
For each file format there is a namedtuple (IOFormat) that has the following
elements:

* a read(filename, index, **kwargs) generator that will yield Atoms objects
* a write(filename, images) function
* a 'single' boolean (False if multiple configurations is supported)
* a 'acceptsfd' boolean (True if file-descriptors are accepted)

There is a dict 'ioformats' that is filled with IOFormat objects as they are
needed.  The 'initialize()' function will create the IOFormat object by
looking at the all_formats dict and by importing the correct read/write
functions from the correct module.  The 'single' and 'acceptsfd' bools are
parsed from two-charcter string in the all_formats dict below.


Example
=======

The xyz format is implemented in the ase/io/xyz.py file which has a
read_xyz() generator and a write_xyz() function.

"""

import functools
import inspect
import os
import sys
from pathlib import Path, PurePath

from ase.atoms import Atoms
from ase.utils import import_module, basestring
from ase.parallel import parallel_function, parallel_generator


class UnknownFileTypeError(Exception):
    pass


class IOFormat:
    def __init__(self, name, desc, code, module_name, encoding=None):
        self.name = name
        self.description = desc
        assert len(code) == 2
        assert code[0] in list('+1')
        assert code[1] in list('BFS')
        self.code = code
        self.module_name = module_name
        self.encoding = encoding

        # (To be set by define_io_format())
        self.extensions = []
        self.globs = []
        self.magic = []

    def open(self, fname, mode='r'):
        # We might want append mode, too
        # We can allow more flags as needed (buffering etc.)
        if mode not in list('rwa'):
            raise ValueError("Only modes allowed are 'r', 'w', and 'a'")
        if mode == 'r' and self.can_read:
            raise NotImplementedError('No reader implemented for {} format'
                                      .format(self.name))
        if mode == 'w' and self.can_write:
            raise NotImplementedError('No writer implemented for {} format'
                                      .format(self.name))
        if mode == 'a' and not self.can_append:
            raise NotImplementedError('Appending not supported by {} format'
                                      .format(self.name))

        if self.isbinary:
            mode += 'b'

        path = Path(fname)
        return path.open(mode, encoding=self.encoding)

    @property
    def can_read(self):
        return self.read is not None

    @property
    def can_write(self):
        return self.write is not None

    @property
    def can_append(self):
        return self.writable and 'append' in self.write.__code__.co_varnames

    def __repr__(self):
        tokens = ['{}={}'.format(name, repr(value))
                  for name, value in vars(self).items()]
        return 'IOFormat({})'.format(', '.join(tokens))

    def __getitem__(self, i):
        return (self.description, self.code)[i]

    @property
    def single(self):
        return self.code[0] == '1'

    @property
    def _formatname(self):
        return self.name.replace('-', '_')

    @property
    def read(self):
        read = getattr(self.module, 'read_' + self._formatname, None)
        if read and not inspect.isgeneratorfunction(read):
            read = functools.partial(wrap_read_function, read)
        return read

    @property
    def write(self):
        return getattr(self.module, 'write_' + self._formatname, None)

    @property
    def modes(self):
        modes = ''
        if self.read:
            modes += 'r'
        if self.write:
            modes += 'w'
        return modes

    def full_description(self):
        lines = ['Name:        {name}',
                 'Description: {description}',
                 'Modes:       {modes}',
                 'Encoding:    {encoding}',
                 'Module:      {module_name}',
                 'Code:        {code}',
                 'Extensions:  {extensions}',
                 'Globs:       {globs}',
                 'Magic:       {magic}']
        desc = '\n'.join(lines)

        myvars = {name: getattr(self, name) for name in dir(self)}
        return desc.format(**myvars)

    @property
    def acceptsfd(self):
        return self.code[1] != 'S'

    @property
    def isbinary(self):
        return self.code[1] == 'B'

    @property
    def module(self):
        if not self.module_name.startswith('ase.io.'):
            raise ValueError('Will only import modules from ase.io, '
                             'not {}'.format(self.module_name))
        try:
            return import_module(self.module_name)
        except ImportError as err:
            raise UnknownFileTypeError('File format not recognized: %s.  '
                                       'Error: %s' % (format, err))

    def match_name(self, basename):
        from fnmatch import fnmatch
        return any(fnmatch(basename, pattern)
                   for pattern in self.globs)

    def match_magic(self, data):
        from fnmatch import fnmatchcase
        return any(fnmatchcase(data, magic + b'*')
                   for magic in self.magic)


ioformats = {}  # will be filled at run-time

# 1=single, +=multiple, F=accepts a file-descriptor, S=needs a file-name str,
# B=like F, but opens in binary mode
all_formats = ioformats  # XXX We should keep one of these.

#glob_patterns = {}
format2modulename = {}  # Left for compatibility only.  Please do not use.
extension2format = {}

def define_io_format(name, desc, code, *, module=None, ext=None,
                     glob=None, magic=None, encoding=None):
    if module is None:
        module = name.replace('-', '_')
        format2modulename[name] = module

    def normalize_patterns(strings):
        if strings is None:
            strings = []
        elif isinstance(strings, (str, bytes)):
            strings = [strings]
        else:
            strings = list(strings)
        return strings

    fmt = IOFormat(name, desc, code, module_name='ase.io.' + module,
                   encoding=encoding)
    fmt.extensions = normalize_patterns(ext)
    fmt.globs = normalize_patterns(glob)
    fmt.magic = normalize_patterns(magic)

    for ext in fmt.extensions:
        if ext in extension2format:
            raise ValueError('extension "{}" already registered'.format(ext))
        extension2format[ext] = fmt

    ioformats[name] = fmt
    return fmt

def get_ioformat(name):
    # This function is left only for backwards compatibility.
    return ioformats[name]

F = define_io_format
F('abinit', 'ABINIT input file', '1F'),
F('aims', 'FHI-aims geometry file', '1S',ext='in'),
F('aims-output', 'FHI-aims output', '+S',
  module='aims', magic=b'*Invoking FHI-aims ...'),
F('bundletrajectory', 'ASE bundle trajectory', '+S'),
F('castep-castep', 'CASTEP output file', '+F',
  module='castep', ext='castep'),
F('castep-cell', 'CASTEP geom file', '1F',
  module='castep', ext='cell'),
F('castep-geom', 'CASTEP trajectory file', '+F',
  module='castep', ext='geom'),
F('castep-md', 'CASTEP molecular dynamics file', '+F',
  module='castep', ext='md'),
F('castep-phonon', 'CASTEP phonon file', '1F',
  module='castep', ext='phonon'),
F('cfg', 'AtomEye configuration', '1F'),
F('cif', 'CIF-file', '+B'),
F('cmdft', 'CMDFT-file', '1F', glob='*I_info'),
F('cp2k-dcd', 'CP2K DCD file', '+B',
  module='cp2k', ext='dcd'),
F('crystal', 'Crystal fort.34 format', '1S',
  ext=['f34', '34'], glob=['f34', '34']),
F('cube', 'CUBE file', '1F'),
F('dacapo', 'Dacapo netCDF output file', '1F'),
F('dacapo-text', 'Dacapo text output', '1F',
  module='dacapo', magic=b'*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n'),
F('db', 'ASE SQLite database file', '+S'),
F('dftb', 'DftbPlus input file', '1S', magic=b'Geometry'),
F('dlp4', 'DL_POLY_4 CONFIG file', '1F',
  module='dlp4', ext='config', glob=['*CONFIG*']),
F('dlp-history', 'DL_POLY HISTORY file', '+F',
  module='dlp4', glob='HISTORY'),
F('dmol-arc', 'DMol3 arc file', '+S',
  module='dmol'),
F('dmol-car', 'DMol3 structure file', '1S',
  module='dmol', ext='car'),
F('dmol-incoor', 'DMol3 structure file', '1S',
  module='dmol'),
F('elk', 'ELK atoms definition', '1S'),
F('eon', 'EON CON file', '+F',
  ext='con'),
F('eps', 'Encapsulated Postscript', '1S'),
F('espresso-in', 'Quantum espresso in file', '1F',
  module='espresso', ext='pwi', magic=[b'*\n&system', b'*\n&SYSTEM']),
F('espresso-out', 'Quantum espresso out file', '+F',
  module='espresso', ext=['out', 'pwo'], magic=b'*Program PWSCF'),
F('etsf', 'ETSF format', '1S'),
F('exciting', 'exciting input', '1S',glob='input.xml'),
F('extxyz', 'Extended XYZ file', '+F'),
F('findsym', 'FINDSYM-format', '+F'),
F('gaussian', 'Gaussian com (input) file', '1S',
  ext=['com', 'gjf']),
F('gaussian-out', 'Gaussian output file', '1F',
  module='gaussian', ext='log'),
F('acemolecule-out', 'ACE output file', '1S',
  module='acemolecule'),
F('acemolecule-input', 'ACE input file', '1S',
  module='acemolecule'),
F('gen', 'DFTBPlus GEN format', '1F'),
F('gif', 'Graphics interchange format', '+S',
  module='animation'),
F('gpaw-out', 'GPAW text output', '+F',
  magic=b'*  ___ ___ ___ _ _ _'),
F('gpw', 'GPAW restart-file', '1S',
  magic=[b'- of UlmGPAW', b'AFFormatGPAW']),
F('gromacs', 'Gromacs coordinates', '1S',
  ext='gro'),
F('gromos', 'Gromos96 geometry file', '1F', ext='g96'),
F('html', 'X3DOM HTML', '1F', module='x3d'),
F('iwm', '?', '1F', glob='atoms.dat'),
F('json', 'ASE JSON database file', '+F', module='db'),
F('jsv', 'JSV file format', '1F'),
F('lammps-dump-text', 'LAMMPS text dump file', '+F',
  module='lammpsrun', magic=b'*\nITEM: TIMESTEP\n'),
F('lammps-dump-binary', 'LAMMPS binary dump file', '+B',
  module='lammpsrun')
F('lammps-data', 'LAMMPS data file', '1F', module='lammpsdata',
  encoding='ascii'
),
F('magres', 'MAGRES ab initio NMR data file', '1F'),
F('mol', 'MDL Molfile', '1F'),
F('mp4', 'MP4 animation', '+S',
  module='animation'),
F('mustem', 'muSTEM xtl file', '1F',
  ext='xtl'),
F('mysql', 'ASE MySQL database file', '+S',
  module='db'),
F('netcdftrajectory', 'AMBER NetCDF trajectory file', '+S'),
F('nomad-json', 'JSON from Nomad archive', '+F',
  ext='nomad-json'),
F('nwchem-in', 'NWChem input file', '1F',
  module='nwchem', ext='nwi'),
F('nwchem-out', 'NWChem output file', '+F',
  module='nwchem', ext='nwo',
  magic=b'*Northwest Computational Chemistry Package'),
F('octopus', 'Octopus input file', '1F', glob='inp'),
F('proteindatabank', 'Protein Data Bank', '+F',
  ext='pdb'),
F('png', 'Portable Network Graphics', '1B'),
F('postgresql', 'ASE PostgreSQL database file', '+S', module='db'),
F('pov', 'Persistance of Vision', '1S'),
F('py', 'Python file', '+F'),
F('qbox', 'QBOX output file', '+F',
  magic=b'*:simulation xmlns:'),
F('res', 'SHELX format', '1S', ext='shelx'),
F('rmc6f', 'RMCProfile', '1S', ext='rmc6f'),
F('sdf', 'SDF format', '1F'),
F('struct', 'WIEN2k structure file', '1S', module='wien2k'),
F('struct_out', 'SIESTA STRUCT file', '1F', module='siesta'),
F('traj', 'ASE trajectory', '+B', module='trajectory',
  magic=[b'- of UlmASE-Trajectory', b'AFFormatASE-Trajectory']),
F('trj', 'Old ASE pickle trajectory', '+S',
  module='pickletrajectory', magic=b'PickleTrajectory'),
F('turbomole', 'TURBOMOLE coord file', '1F', glob='coord',
  magic=b'$coord'),
F('turbomole-gradient', 'TURBOMOLE gradient file', '+F',
  module='turbomole', glob='gradient', magic=b'$grad'),
F('v-sim', 'V_Sim ascii file', '1F', ext='ascii'),
F('vasp', 'VASP POSCAR/CONTCAR', '1F',
  ext='poscar', glob=['*POSCAR*', '*CONTCAR*']),
F('vasp-out', 'VASP OUTCAR file', '+F',
  module='vasp', glob='*OUTCAR*'),
F('vasp-xdatcar', 'VASP XDATCAR file', '+F',
  module='vasp', glob='*XDATCAR*'),
F('vasp-xml', 'VASP vasprun.xml file', '+F',
  module='vasp', glob='*vasp*.xml'),
F('vti', 'VTK XML Image Data', '1F', module='vtkxml'),
F('vtu', 'VTK XML Unstructured Grid', '1F', module='vtkxml'),
F('x3d', 'X3D', '1S'),
F('xsd', 'Materials Studio file', '1F'),
F('xsf', 'XCrySDen Structure File', '+F',
  magic=[b'*\nANIMSTEPS', b'*\nCRYSTAL', b'*\nSLAB', b'*\nPOLYMER',
         b'*\nMOLECULE', b'*\nATOMS']),
F('xtd', 'Materials Studio file', '+F'),
F('xyz', 'XYZ-file', '+F')

netcdfconventions2format = {
    'http://www.etsf.eu/fileformats': 'etsf',
    'AMBER': 'netcdftrajectory'
}


def get_compression(filename):
    """
    Parse any expected file compression from the extension of a filename.
    Return the filename without the extension, and the extension. Recognises
    ``.gz``, ``.bz2``, ``.xz``.

    >>> get_compression('H2O.pdb.gz')
    ('H2O.pdb', 'gz')
    >>> get_compression('crystal.cif')
    ('crystal.cif', None)

    Parameters
    ==========
    filename: str
        Full filename including extension.

    Returns
    =======
    (root, extension): (str, str or None)
        Filename split into root without extension, and the extension
        indicating compression format. Will not split if compression
        is not recognised.
    """
    # Update if anything is added
    valid_compression = ['gz', 'bz2', 'xz']

    # Use stdlib as it handles most edge cases
    root, compression = os.path.splitext(filename)

    # extension keeps the '.' so remember to remove it
    if compression.strip('.') in valid_compression:
        return root, compression.strip('.')
    else:
        return filename, None


def open_with_compression(filename, mode='r'):
    """
    Wrapper around builtin `open` that will guess compression of a file
    from the filename and open it for reading or writing as if it were
    a standard file.

    Implemented for ``gz``(gzip), ``bz2``(bzip2) and ``xz``(lzma).

    Supported modes are:
       * 'r', 'rt', 'w', 'wt' for text mode read and write.
       * 'rb, 'wb' for binary read and write.

    Parameters
    ==========
    filename: str
        Path to the file to open, including any extensions that indicate
        the compression used.
    mode: str
        Mode to open the file, same as for builtin ``open``, e.g 'r', 'w'.

    Returns
    =======
    fd: file
        File-like object open with the specified mode.
    """

    # Compressed formats sometimes default to binary, so force text mode.
    if mode == 'r':
        mode = 'rt'
    elif mode == 'w':
        mode = 'wt'
    elif mode == 'a':
        mode = 'at'

    root, compression = get_compression(filename)

    if compression is None:
        return open(filename, mode)
    elif compression == 'gz':
        import gzip
        fd = gzip.open(filename, mode=mode)
    elif compression == 'bz2':
        import bz2
        fd = bz2.open(filename, mode=mode)
    elif compression == 'xz':
        try:
            from lzma import open as lzma_open
        except ImportError:
            from backports.lzma import open as lzma_open
        fd = lzma_open(filename, mode)
    else:
        fd = open(filename, mode)

    return fd


def wrap_read_function(read, filename, index=None, **kwargs):
    """Convert read-function to generator."""
    if index is None:
        yield read(filename, **kwargs)
    else:
        for atoms in read(filename, index, **kwargs):
            yield atoms


def write(filename, images, format=None, parallel=True, append=False,
          **kwargs):
    """Write Atoms object(s) to file.

    filename: str or file
        Name of the file to write to or a file descriptor.  The name '-'
        means standard output.
    images: Atoms object or list of Atoms objects
        A single Atoms object or a list of Atoms objects.
    format: str
        Used to specify the file-format.  If not given, the
        file-format will be taken from suffix of the filename.
    parallel: bool
        Default is to write on master only.  Use parallel=False to write
        from all slaves.
    append: bool
        Default is to open files in 'w' or 'wb' mode, overwriting
        existing files.  In some cases opening the file in 'a' or 'ab'
        mode (appending) is usefull,
        e.g. writing trajectories or saving multiple Atoms objects in one file.
        WARNING: If the file format does not support multiple entries without
        additional keywords/headers, files created using 'append=True'
        might not be readable by any program! They will nevertheless be
        written without error message.

    The use of additional keywords is format specific."""

    if isinstance(filename, PurePath):
        filename = str(filename)

    if isinstance(filename, basestring):
        filename = os.path.expanduser(filename)
        fd = None
        if filename == '-':
            fd = sys.stdout
            filename = None
        elif format is None:
            format = filetype(filename, read=False)
            assert isinstance(format, str)
    else:
        fd = filename
        filename = None

    format = format or 'json'  # default is json

    io = ioformats[format]

    _write(filename, fd, format, io, images, parallel=parallel, append=append,
           **kwargs)


@parallel_function
def _write(filename, fd, format, io, images, parallel=None, append=False,
           **kwargs):
    if isinstance(images, Atoms):
        images = [images]

    if io.single:
        if len(images) > 1:
            raise ValueError('{}-format can only store 1 Atoms object.'
                             .format(format))
        images = images[0]

    if io.write is None:
        raise ValueError("Can't write to {}-format".format(format))

    # Special case for json-format:
    if format == 'json' and (len(images) > 1 or append):
        if filename is not None:
            io.write(filename, images, append=append, **kwargs)
            return
        raise ValueError("Can't write more than one image to file-descriptor "
                         'using json-format.')

    if io.acceptsfd:
        open_new = (fd is None)
        if open_new:
            mode = 'wb' if io.isbinary else 'w'
            if append:
                mode = mode.replace('w', 'a')
            fd = open_with_compression(filename, mode)
            # XXX remember to re-enable compressed open
            #fd = io.open(filename, mode)
        io.write(fd, images, **kwargs)
        if open_new:
            fd.close()
    else:
        if fd is not None:
            raise ValueError("Can't write {}-format to file-descriptor"
                             .format(format))
        if 'append' in io.write.__code__.co_varnames:
            io.write(filename, images, append=append, **kwargs)
        elif append:
            raise ValueError("Cannot append to {}-format, write-function "
                             "does not support the append keyword."
                             .format(format))
        else:
            io.write(filename, images, **kwargs)


def read(filename, index=None, format=None, parallel=True, **kwargs):
    """Read Atoms object(s) from file.

    filename: str or file
        Name of the file to read from or a file descriptor.
    index: int, slice or str
        The last configuration will be returned by default.  Examples:

            * ``index=0``: first configuration
            * ``index=-2``: second to last
            * ``index=':'`` or ``index=slice(None)``: all
            * ``index='-3:'`` or ``index=slice(-3, None)``: three last
            * ``index='::2'`` or ``index=slice(0, None, 2)``: even
            * ``index='1::2'`` or ``index=slice(1, None, 2)``: odd
    format: str
        Used to specify the file-format.  If not given, the
        file-format will be guessed by the *filetype* function.
    parallel: bool
        Default is to read on master and broadcast to slaves.  Use
        parallel=False to read on all slaves.

    Many formats allow on open file-like object to be passed instead
    of ``filename``. In this case the format cannot be auto-decected,
    so the ``format`` argument should be explicitly given."""

    if isinstance(filename, PurePath):
        filename = str(filename)
    if filename == '-':
        filename = sys.stdin
    if isinstance(index, basestring):
        try:
            index = string2index(index)
        except ValueError:
            pass

    filename, index = parse_filename(filename, index)
    if index is None:
        index = -1
    format = format or filetype(filename)

    io = ioformats[format]
    if isinstance(index, (slice, basestring)):
        return list(_iread(filename, index, format, io, parallel=parallel,
                           **kwargs))
    else:
        return next(_iread(filename, slice(index, None), format, io,
                           parallel=parallel, **kwargs))


def iread(filename, index=None, format=None, parallel=True, **kwargs):
    """Iterator for reading Atoms objects from file.

    Works as the `read` function, but yields one Atoms object at a time
    instead of all at once."""

    if isinstance(index, basestring):
        index = string2index(index)

    filename, index = parse_filename(filename, index)

    if index is None or index == ':':
        index = slice(None, None, None)

    if not isinstance(index, (slice, basestring)):
        index = slice(index, (index + 1) or None)

    format = format or filetype(filename)
    io = ioformats[format]

    for atoms in _iread(filename, index, format, io, parallel=parallel,
                        **kwargs):
        yield atoms


@parallel_generator
def _iread(filename, index, format, io, parallel=None, full_output=False,
           **kwargs):
    if isinstance(filename, basestring):
        filename = os.path.expanduser(filename)

    if not io.read:
        raise ValueError("Can't read from {}-format".format(format))

    if io.single:
        start = index.start
        assert start is None or start == 0 or start == -1
        args = ()
    else:
        args = (index,)

    must_close_fd = False
    if isinstance(filename, basestring):
        if io.acceptsfd:
            mode = 'rb' if io.isbinary else 'r'
            fd = open_with_compression(filename, mode)
            must_close_fd = True
        else:
            fd = filename
    else:
        assert io.acceptsfd
        fd = filename

    # Make sure fd is closed in case loop doesn't finish:
    try:
        for dct in io.read(fd, *args, **kwargs):
            if not isinstance(dct, dict):
                dct = {'atoms': dct}
            if full_output:
                yield dct
            else:
                yield dct['atoms']
    finally:
        if must_close_fd:
            fd.close()


def parse_filename(filename, index=None):
    if not isinstance(filename, basestring):
        return filename, index

    extension = os.path.basename(filename)
    if '@' not in extension:
        return filename, index

    newindex = None
    newfilename, newindex = filename.rsplit('@', 1)

    if isinstance(index, slice):
        return newfilename, index
    try:
        newindex = string2index(newindex)
    except ValueError:
        pass

    return newfilename, newindex


def string2index(string):
    """Convert index string to either int or slice"""
    if ':' not in string:
        return int(string)
    i = []
    for s in string.split(':'):
        if s == '':
            i.append(None)
        else:
            i.append(int(s))
    i += (3 - len(i)) * [None]
    return slice(*i)


def filetype(filename, read=True, guess=True):
    """Try to guess the type of the file.

    First, special signatures in the filename will be checked for.  If that
    does not identify the file type, then the first 2000 bytes of the file
    will be read and analysed.  Turn off this second part by using
    read=False.

    Can be used from the command-line also::

        $ ase info filename ...
    """

    ext = None
    if isinstance(filename, basestring):
        if os.path.isdir(filename):
            if os.path.basename(os.path.normpath(filename)) == 'states':
                return 'eon'
            return 'bundletrajectory'

        if filename.startswith('postgres'):
            return 'postgresql'

        if filename.startswith('mysql') or filename.startswith('mariadb'):
            return 'mysql'

        # strip any compression extensions that can be read
        root, compression = get_compression(filename)
        basename = os.path.basename(root)

        if '.' in basename:
            ext = os.path.splitext(basename)[1].strip('.').lower()
            if ext in ['xyz', 'cube', 'json', 'cif']:
                return ext

        for fmt in ioformats.values():
            if fmt.match_name(basename):
                return fmt.name


        if not read:
            if ext is None:
                raise UnknownFileTypeError('Could not guess file type')
            ioformat = extension2format.get(ext)
            if ioformat:
                return ioformat.name

            # askhl: This is strange, we don't know if ext is a format:
            return ext

        fd = open_with_compression(filename, 'rb')
    else:
        fd = filename
        if fd is sys.stdin:
            return 'json'

    data = fd.read(50000)
    if fd is not filename:
        fd.close()
    else:
        fd.seek(0)

    if len(data) == 0:
        raise UnknownFileTypeError('Empty file: ' + filename)

    if data.startswith(b'CDF'):
        # We can only recognize these if we actually have the netCDF4 module.
        try:
            import netCDF4
        except ImportError:
            pass
        else:
            nc = netCDF4.Dataset(filename)
            if 'Conventions' in nc.ncattrs():
                if nc.Conventions in netcdfconventions2format:
                    return netcdfconventions2format[nc.Conventions]
                else:
                    raise UnknownFileTypeError(
                        "Unsupported NetCDF convention: "
                        "'{}'".format(nc.Conventions))
            else:
                raise UnknownFileTypeError("NetCDF file does not have a "
                                           "'Conventions' attribute.")

    for ioformat in ioformats.values():
        if ioformat.match_magic(data):
            return ioformat.name

    format = None
    if ext in extension2format:
        format = extension2format[ext].name

    if format is None and guess:
        format = ext
    if format is None:
        # Do quick xyz check:
        lines = data.splitlines()
        if lines and lines[0].strip().isdigit():
            return 'xyz'

        raise UnknownFileTypeError('Could not guess file type')
    assert isinstance(format, str)
    return format


def index2range(index, nsteps):
    """Method to convert a user given *index* option to a list of indices.

    Returns a range.
    """
    if isinstance(index, int):
        if index < 0:
            tmpsnp = nsteps + index
            trbl = range(tmpsnp, tmpsnp + 1, 1)
        else:
            trbl = range(index, index + 1, 1)
    elif isinstance(index, slice):
        start = index.start
        stop = index.stop
        step = index.step

        if start is None:
            start = 0
        elif start < 0:
            start = nsteps + start

        if step is None:
            step = 1

        if stop is None:
            stop = nsteps
        elif stop < 0:
            stop = nsteps + stop

        trbl = range(start, stop, step)
    else:
        raise RuntimeError("index2range handles integers and slices only.")
    return trbl
