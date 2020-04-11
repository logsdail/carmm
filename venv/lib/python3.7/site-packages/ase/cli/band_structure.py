
import json
import numpy as np

from ase.io import read
from ase.io.formats import UnknownFileTypeError
from ase.io.jsonio import read_json
from ase.geometry import crystal_structure_from_cell
from ase.dft.kpoints import (get_monkhorst_pack_size_and_offset,
                             monkhorst_pack_interpolate,
                             bandpath, BandPath,
                             get_special_points)
from ase.dft.band_structure import BandStructure


class CLICommand:
    """Plot band-structure.

    Read eigenvalues and k-points from file and plot result from
    band-structure calculation or interpolate
    from Monkhorst-Pack sampling to a given path (--path=PATH).

    Example:

        $ ase band-structure al.gpw -r -10 10
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('calculation',
                            help='Path to output file(s) from calculation.')
        parser.add_argument('-q', '--quiet', action='store_true',
                            help='Less output.')
        parser.add_argument('-k', '--path', help='Example "GXL".')
        parser.add_argument('-n', '--points', type=int, default=100,
                            help='Number of points along the path '
                            '(default: 100)')
        parser.add_argument('-r', '--range', nargs=2, default=['-3', '3'],
                            metavar=('emin', 'emax'),
                            help='Default: "-3.0 3.0" '
                            '(in eV relative to Fermi level).')

    @staticmethod
    def run(args, parser):
        main(args, parser)

def atoms2bandstructure(atoms, parser, args):
    cell = atoms.get_cell()
    calc = atoms.calc
    bzkpts = calc.get_bz_k_points()
    ibzkpts = calc.get_ibz_k_points()
    efermi = calc.get_fermi_level()
    nibz = len(ibzkpts)
    nspins = 1 + int(calc.get_spin_polarized())

    eps = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                     for k in range(nibz)]
                    for s in range(nspins)])
    if not args.quiet:
        print('Spins, k-points, bands: {}, {}, {}'.format(*eps.shape))

    if bzkpts is None:
        if ibzkpts is None:
            raise ValueError('Cannot find any k-point data')
        else:
            path_kpts = ibzkpts
    else:
        try:
            size, offset = get_monkhorst_pack_size_and_offset(bzkpts)
        except ValueError:
            path_kpts = ibzkpts
        else:
            if not args.quiet:
                print('Interpolating from Monkhorst-Pack grid (size, offset):')
                print(size, offset)
            if args.path is None:
                err = 'Please specify a path!'
                try:
                    cs = crystal_structure_from_cell(cell)
                except ValueError:
                    err += ('\nASE cannot automatically '
                            'recognize this crystal structure')
                else:
                    from ase.dft.kpoints import special_paths
                    kptpath = special_paths[cs]
                    err += ('\nIt looks like you have a {} crystal structure.'
                            '\nMaybe you want its special path:'
                            ' {}'.format(cs, kptpath))
                parser.error(err)
        bz2ibz = calc.get_bz_to_ibz_map()

        path_kpts = bandpath(args.path, atoms.cell, args.points)[0]
        icell = atoms.get_reciprocal_cell()
        eps = monkhorst_pack_interpolate(path_kpts, eps.transpose(1, 0, 2),
                                         icell, bz2ibz, size, offset)
        eps = eps.transpose(1, 0, 2)

    special_points = get_special_points(cell)
    path = BandPath(atoms.cell, kpts=path_kpts,
                    special_points=special_points)

    return BandStructure(path, eps, reference=efermi)


def read_band_structure(args, parser):
    # Read first as Atoms object, then try as JSON band structure:
    try:
        atoms = read(args.calculation)
    except UnknownFileTypeError:
        pass
    else:
        return atoms2bandstructure(atoms, parser, args)

    try:
        bs = read_json(args.calculation)
    except json.decoder.JSONDecodeError:
        parser.error('File resembles neither atoms nor band structure')

    objtype = getattr(bs, 'ase_objtype', None)
    if objtype != 'bandstructure':
        parser.error('Expected band structure, but this file contains a {}'
                     .format(objtype or type(bs)))

    return bs

def main(args, parser):
    import matplotlib.pyplot as plt
    bs = read_band_structure(args, parser)
    emin, emax = (float(e) for e in args.range)
    fig = plt.gcf()
    fig.canvas.set_window_title(args.calculation)
    ax = fig.gca()
    bs.plot(ax=ax, emin=emin, emax=emax)
    plt.show()
