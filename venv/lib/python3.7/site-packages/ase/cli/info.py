import platform
import os
import sys

from ase.utils import import_module
from ase.utils import search_current_git_hash
from ase.io.formats import filetype, ioformats, UnknownFileTypeError
from ase.io.ulm import print_ulm_info
from ase.io.bundletrajectory import print_bundletrajectory_info


class CLICommand:
    """Print information about files or system.

    Without any filename(s), informations about the ASE installation will be
    shown (Python version, library versions, ...).

    With filename(s), the file format will be determined for each file.
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('filename', nargs='*',
                            help='Name of file to determine format for.')
        parser.add_argument('-v', '--verbose', action='store_true',
                            help='Show more information about files.')
        parser.add_argument('--formats', action='store_true',
                            help='List file formats known to ASE.')
        parser.add_argument('--calculators', action='store_true',
                            help='List calculators known to ASE '
                            'and whether they appear to be installed.')

    @staticmethod
    def run(args):
        if not args.filename:
            print_info()
            if args.formats:
                print()
                print_formats()
            if args.calculators:
                print()
                from ase.calculators.autodetect import (detect_calculators,
                                                        format_configs)
                configs = detect_calculators()
                print('Calculators:')
                for message in format_configs(configs):
                    print('  {}'.format(message))
                print()
                print('Available: {}'.format(','.join(sorted(configs))))
            return

        n = max(len(filename) for filename in args.filename) + 2
        for filename in args.filename:
            try:
                format = filetype(filename)
            except FileNotFoundError:
                format = '?'
                description = 'No such file'
            except UnknownFileTypeError:
                format = '?'
                description = '?'
            else:
                if format in ioformats:
                    description = ioformats[format].description
                else:
                    description = '?'

            print('{:{}}{} ({})'.format(filename + ':', n,
                                        description, format))
            if args.verbose:
                if format == 'traj':
                    print_ulm_info(filename)
                elif format == 'bundletrajectory':
                    print_bundletrajectory_info(filename)


def print_info():
    versions = [('platform', platform.platform()),
                ('python-' + sys.version.split()[0], sys.executable)]
    for name in ['ase', 'numpy', 'scipy', 'ase_ext']:
        try:
            module = import_module(name)
        except ImportError:
            if name != 'ase_ext':
                versions.append((name, 'no'))
        else:
            # Search for git hash
            githash = search_current_git_hash(module)
            if githash is None:
                githash = ''
            else:
                githash = '-{:.10}'.format(githash)
            versions.append((name + '-' + module.__version__ + githash,
                            module.__file__.rsplit(os.sep, 1)[0] + os.sep))

    for a, b in versions:
        print('{:25}{}'.format(a, b))


def print_formats():
    print('Supported formats:')
    for fmtname in sorted(ioformats):
        fmt = ioformats[fmtname]

        infos = [fmt.modes, 'single' if fmt.single else 'multi']
        if fmt.isbinary:
            infos.append('binary')
        if fmt.encoding is not None:
            infos.append(fmt.encoding)
        infostring = '/'.join(infos)

        moreinfo = [infostring]
        if fmt.extensions:
            moreinfo.append('ext={}'.format('|'.join(fmt.extensions)))
        if fmt.globs:
            moreinfo.append('glob={}'.format('|'.join(fmt.globs)))

        print('  {} [{}]: {}'.format(fmt.name,
                                     ', '.join(moreinfo),
                                     fmt.description))
