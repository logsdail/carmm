"""
This module defines the ASE interface to SIESTA.

Written by Mads Engelund (see www.espeem.com)

Home of the SIESTA package:
http://www.uam.es/departamentos/ciencias/fismateriac/siesta

2017.04 - Pedro Brandimarte: changes for python 2-3 compatible

"""

import os
import re
import tempfile
import warnings
import shutil
from os.path import join, isfile, islink
import numpy as np

from ase.units import Ry, eV, Bohr
from ase.data import atomic_numbers
from ase.calculators.siesta.import_functions import read_rho, xv_to_atoms
from ase.calculators.siesta.import_functions import \
    get_valence_charge, read_vca_synth_block
from ase.calculators.calculator import FileIOCalculator, ReadError
from ase.calculators.calculator import Parameters, all_changes
from ase.calculators.siesta.parameters import PAOBasisBlock, Species
from ase.calculators.siesta.parameters import format_fdf


meV = 0.001 * eV


def get_siesta_version(command):
    """ Return SIESTA version number.

    Run the command, for instance 'siesta' and
    then parse the output in order find the
    version number.
    """

    temp_dirname = tempfile.mkdtemp(prefix='siesta-version-check-')
    try:
        from subprocess import Popen, PIPE
        proc = Popen(command,
                     stdin=PIPE,
                     stdout=PIPE,
                     stderr=PIPE,
                     cwd=temp_dirname)
        output, _ = proc.communicate()
        # We are not providing any input, so Siesta will give us a failure
        # saying that it has no Chemical_species_label and exit status 1
        # (as of siesta-4.1-b4)
    finally:
        shutil.rmtree(temp_dirname)

    match = re.search(rb'Siesta Version\s*:\s*(\S+)', output)

    if match is None:
        raise RuntimeError('Could not get Siesta version info from output '
                           '{!r}'.format(output))

    string = match.group(1).decode('ascii')
    return string


def bandpath2bandpoints(path):
    lines = []
    add = lines.append

    add('BandLinesScale ReciprocalLatticeVectors\n')
    add('%block BandPoints\n')
    for kpt in path.kpts:
        add('    {:18.15f} {:18.15f} {:18.15f}\n'.format(*kpt))
    add('%endblock BandPoints')
    return ''.join(lines)


def read_bands_file(fd):
    efermi = float(next(fd))
    next(fd)  # Appears to be max/min energy.  Not important for us
    header = next(fd)  # Array shape: nbands, nspins, nkpoints
    nbands, nspins, nkpts = np.array(header.split()).astype(int)

    # three fields for kpt coords, then all the energies
    ntokens = nbands * nspins + 3

    # Read energies for each kpoint:
    data = []
    for i in range(nkpts):
        line = next(fd)
        tokens = line.split()
        while len(tokens) < ntokens:
            # Multirow table.  Keep adding lines until the table ends,
            # which should happen exactly when we have all the energies
            # for this kpoint.
            line = next(fd)
            tokens += line.split()
        assert len(tokens) == ntokens
        values = np.array(tokens).astype(float)
        data.append(values)

    data = np.array(data)
    assert len(data) == nkpts
    kpts = data[:, :3]
    energies = data[:, 3:]
    energies = energies.reshape(nkpts, nspins, nbands)
    assert energies.shape == (nkpts, nspins, nbands)
    return kpts, energies, efermi


def resolve_band_structure(path, kpts, energies, efermi):
    """Convert input BandPath along with Siesta outputs into BS object."""
    # Right now this function doesn't do much.
    #
    # Not sure how the output kpoints in the siesta.bands file are derived.
    # They appear to be related to the lattice parameter.
    #
    # We should verify that they are consistent with our input path,
    # but since their meaning is unclear, we can't quite do so.
    #
    # Also we should perhaps verify the cell.  If we had the cell, we
    # could construct the bandpath from scratch (i.e., pure outputs).
    from ase.dft.band_structure import BandStructure
    ksn2e = energies
    skn2e = np.swapaxes(ksn2e, 0, 1)
    bs = BandStructure(path, skn2e, reference=efermi)
    return bs


class SiestaParameters(Parameters):
    """Parameters class for the calculator.
    Documented in BaseSiesta.__init__

    """
    def __init__(
            self,
            label='siesta',
            mesh_cutoff=200 * Ry,
            energy_shift=100 * meV,
            kpts=None,
            xc='LDA',
            basis_set='DZP',
            spin='non-polarized',
            species=tuple(),
            pseudo_qualifier=None,
            pseudo_path=None,
            symlink_pseudos=None,
            atoms=None,
            restart=None,
            ignore_bad_restart_file=False,
            fdf_arguments=None,
            atomic_coord_format='xyz',
            bandpath=None):
        kwargs = locals()
        kwargs.pop('self')
        Parameters.__init__(self, **kwargs)


class Siesta(FileIOCalculator):
    """Calculator interface to the SIESTA code.
    """
    allowed_basis_names = ['SZ', 'SZP', 'DZ', 'DZP', 'TZP']
    allowed_spins = ['non-polarized', 'collinear',
                     'non-collinear', 'spin-orbit']
    allowed_xc = {
        'LDA': ['PZ', 'CA', 'PW92'],
        'GGA': ['PW91', 'PBE', 'revPBE', 'RPBE',
                'WC', 'AM05', 'PBEsol', 'PBEJsJrLO',
                'PBEGcGxLO', 'PBEGcGxHEG', 'BLYP'],
        'VDW': ['DRSLL', 'LMKLL', 'KBM', 'C09', 'BH', 'VV']}

    name = 'siesta'
    command = 'siesta < PREFIX.fdf > PREFIX.out'
    implemented_properties = (
        'energy',
        'forces',
        'stress',
        'dipole',
        'eigenvalues',
        'density',
        'fermi_energy')

    # Dictionary of valid input vaiables.
    default_parameters = SiestaParameters()

    # XXX Not a ASE standard mechanism (yet).  We need to communicate to
    # ase.dft.band_structure.calculate_band_structure() that we expect
    # it to use the bandpath keyword.
    accepts_bandpath_keyword = True

    def __init__(self, command=None, **kwargs):
        """ASE interface to the SIESTA code.

        Parameters:
           - label        : The basename of all files created during
                            calculation.
           - mesh_cutoff  : Energy in eV.
                            The mesh cutoff energy for determining number of
                            grid points in the matrix-element calculation.
           - energy_shift : Energy in eV
                            The confining energy of the basis set generation.
           - kpts         : Tuple of 3 integers, the k-points in different
                            directions.
           - xc           : The exchange-correlation potential. Can be set to
                            any allowed value for either the Siesta
                            XC.funtional or XC.authors keyword. Default "LDA"
           - basis_set    : "SZ"|"SZP"|"DZ"|"DZP"|"TZP", strings which specify
                            the type of functions basis set.
           - spin         : "non-polarized"|"collinear"|
                            "non-collinear|spin-orbit".
                            The level of spin description to be used.
           - species      : None|list of Species objects. The species objects
                            can be used to to specify the basis set,
                            pseudopotential and whether the species is ghost.
                            The tag on the atoms object and the element is used
                            together to identify the species.
           - pseudo_path  : None|path. This path is where
                            pseudopotentials are taken from.
                            If None is given, then then the path given
                            in $SIESTA_PP_PATH will be used.
           - pseudo_qualifier: None|string. This string will be added to the
                            pseudopotential path that will be retrieved.
                            For hydrogen with qualifier "abc" the
                            pseudopotential "H.abc.psf" will be retrieved.
           - symlink_pseudos: None|bool
                            If true, symlink pseudopotentials
                            into the calculation directory, else copy them.
                            Defaults to true on Unix and false on Windows.
           - atoms        : The Atoms object.
           - restart      : str.  Prefix for restart file.
                            May contain a directory.
                            Default is  None, don't restart.
           - ignore_bad_restart_file: bool.
                            Ignore broken or missing restart file.
                            By default, it is an error if the restart
                            file is missing or broken.
           - fdf_arguments: Explicitly given fdf arguments. Dictonary using
                            Siesta keywords as given in the manual. List values
                            are written as fdf blocks with each element on a
                            separate line, while tuples will write each element
                            in a single line.  ASE units are assumed in the
                            input.
           - atomic_coord_format: "xyz"|"zmatrix", strings to switch between
                            the default way of entering the system's geometry
                            (via the block AtomicCoordinatesAndAtomicSpecies)
                            and a recent method via the block Zmatrix. The
                            block Zmatrix allows to specify basic geometry
                            constrains such as realized through the ASE classes
                            FixAtom, FixedLine and FixedPlane.
        """

        # Put in the default arguments.
        parameters = self.default_parameters.__class__(**kwargs)

        # Call the base class.
        FileIOCalculator.__init__(
            self,
            command=command,
            **parameters)

        # For compatibility with old variable name:
        commandvar = os.environ.get('SIESTA_COMMAND')
        if commandvar is not None:
            warnings.warn('Please use $ASE_SIESTA_COMMAND and not '
                          '$SIESTA_COMMAND, which will be ignored '
                          'in the future.  The new command format will not '
                          'work with the "<%s > %s" specification.  Use '
                          'instead e.g. "ASE_SIESTA_COMMAND=siesta'
                          ' < PREFIX.fdf > PREFIX.out", where PREFIX will '
                          'automatically be replaced by calculator label',
                          np.VisibleDeprecationWarning)
            runfile = self.prefix + '.fdf'
            outfile = self.prefix + '.out'
            try:
                self.command = commandvar % (runfile, outfile)
            except TypeError:
                raise ValueError(
                    "The 'SIESTA_COMMAND' environment must " +
                    "be a format string" +
                    " with two string arguments.\n" +
                    "Example : 'siesta < %s > %s'.\n" +
                    "Got '%s'" % commandvar)

    def __getitem__(self, key):
        """Convenience method to retrieve a parameter as
        calculator[key] rather than calculator.parameters[key]

            Parameters:
                -key       : str, the name of the parameters to get.
        """
        return self.parameters[key]

    def species(self, atoms):
        """Find all relevant species depending on the atoms object and
        species input.

            Parameters :
                - atoms : An Atoms object.
        """
        # For each element use default species from the species input, or set
        # up a default species  from the general default parameters.
        symbols = np.array(atoms.get_chemical_symbols())
        tags = atoms.get_tags()
        species = list(self['species'])
        default_species = [
            s for s in species
            if (s['tag'] is None) and s['symbol'] in symbols]
        default_symbols = [s['symbol'] for s in default_species]
        for symbol in symbols:
            if symbol not in default_symbols:
                spec = Species(symbol=symbol,
                               basis_set=self['basis_set'],
                               tag=None)
                default_species.append(spec)
                default_symbols.append(symbol)
        assert len(default_species) == len(np.unique(symbols))

        # Set default species as the first species.
        species_numbers = np.zeros(len(atoms), int)
        i = 1
        for spec in default_species:
            mask = symbols == spec['symbol']
            species_numbers[mask] = i
            i += 1

        # Set up the non-default species.
        non_default_species = [s for s in species if not s['tag'] is None]
        for spec in non_default_species:
            mask1 = (tags == spec['tag'])
            mask2 = (symbols == spec['symbol'])
            mask = np.logical_and(mask1, mask2)
            if sum(mask) > 0:
                species_numbers[mask] = i
                i += 1
        all_species = default_species + non_default_species

        return all_species, species_numbers

    def set(self, **kwargs):
        """Set all parameters.

            Parameters:
                -kwargs  : Dictionary containing the keywords defined in
                           SiestaParameters.
        """

        # Find not allowed keys.
        default_keys = list(self.__class__.default_parameters)
        offending_keys = set(kwargs) - set(default_keys)
        if len(offending_keys) > 0:
            mess = "'set' does not take the keywords: %s "
            raise ValueError(mess % list(offending_keys))

        # Use the default parameters.
        parameters = self.__class__.default_parameters.copy()
        parameters.update(kwargs)
        kwargs = parameters

        # Check energy inputs.
        for arg in ['mesh_cutoff', 'energy_shift']:
            value = kwargs.get(arg)
            if value is None:
                continue
            if not (isinstance(value, (float, int)) and value > 0):
                mess = "'%s' must be a positive number(in eV), \
                    got '%s'" % (arg, value)
                raise ValueError(mess)

        # Check the basis set input.
        if 'basis_set' in kwargs:
            basis_set = kwargs['basis_set']
            allowed = self.allowed_basis_names
            if not (isinstance(basis_set, PAOBasisBlock) or
                    basis_set in allowed):
                mess = "Basis must be either %s, got %s" % (allowed, basis_set)
                raise ValueError(mess)

        # Check the spin input.
        if 'spin' in kwargs:
            if kwargs['spin'] == 'UNPOLARIZED':
                warnings.warn("The keyword 'UNPOLARIZED' is deprecated,"
                              "and replaced by 'non-polarized'",
                              np.VisibleDeprecationWarning)
                kwargs['spin'] = 'non-polarized'

            spin = kwargs['spin']
            if spin is not None and (spin.lower() not in self.allowed_spins):
                mess = "Spin must be %s, got '%s'" % (self.allowed_spins, spin)
                raise ValueError(mess)

        # Check the functional input.
        xc = kwargs.get('xc', 'LDA')
        if isinstance(xc, (tuple, list)) and len(xc) == 2:
            functional, authors = xc
            if functional.lower() not in [k.lower() for k in self.allowed_xc]:
                mess = "Unrecognized functional keyword: '%s'" % functional
                raise ValueError(mess)

            lsauthorslower = [a.lower() for a in self.allowed_xc[functional]]
            if authors.lower() not in lsauthorslower:
                mess = "Unrecognized authors keyword for %s: '%s'"
                raise ValueError(mess % (functional, authors))

        elif xc in self.allowed_xc:
            functional = xc
            authors = self.allowed_xc[xc][0]
        else:
            found = False
            for key, value in self.allowed_xc.items():
                if xc in value:
                    found = True
                    functional = key
                    authors = xc
                    break

            if not found:
                raise ValueError("Unrecognized 'xc' keyword: '%s'" % xc)
        kwargs['xc'] = (functional, authors)

        # Check fdf_arguments.
        if kwargs['fdf_arguments'] is None:
            kwargs['fdf_arguments'] = {}

        if not isinstance(kwargs['fdf_arguments'], dict):
            raise TypeError("fdf_arguments must be a dictionary.")

        # Call baseclass.
        FileIOCalculator.set(self, **kwargs)

    def set_fdf_arguments(self, fdf_arguments):
        """ Set the fdf_arguments after the initialization of the
            calculator.
        """
        self.validate_fdf_arguments(fdf_arguments)
        FileIOCalculator.set(self, fdf_arguments=fdf_arguments)

    def validate_fdf_arguments(self, fdf_arguments):
        """ Raises error if the fdf_argument input is not a
            dictionary of allowed keys.
        """
        # None is valid
        if fdf_arguments is None:
            return

        # Type checking.
        if not isinstance(fdf_arguments, dict):
            raise TypeError("fdf_arguments must be a dictionary.")

    def calculate(self,
                  atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        """Capture the RuntimeError from FileIOCalculator.calculate
        and add a little debug information from the Siesta output.

        See base FileIocalculator for documentation.
        """

        FileIOCalculator.calculate(
            self,
            atoms=atoms,
            properties=properties,
            system_changes=system_changes)

        # The below snippet would run if calculate() failed but I have
        # disabled it for now since it looks to be just for debugging.
        # --askhl
        """
        # Here a test to check if the potential are in the right place!!!
        except RuntimeError as e:
            try:
                fname = os.path.join(self.directory, self.label+'.out')
                with open(fname, 'r') as f:
                    lines = f.readlines()
                debug_lines = 10
                print('##### %d last lines of the Siesta output' % debug_lines)
                for line in lines[-20:]:
                    print(line.strip())
                print('##### end of siesta output')
                raise e
            except:
                raise e
        """

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input (fdf)-file.
        See calculator.py for further details.

        Parameters:
            - atoms        : The Atoms object to write.
            - properties   : The properties which should be calculated.
            - system_changes : List of properties changed since last run.
        """
        # Call base calculator.
        FileIOCalculator.write_input(
            self,
            atoms=atoms,
            properties=properties,
            system_changes=system_changes)

        if system_changes is None and properties is None:
            return

        filename = self.getpath(ext='fdf')

        # On any changes, remove all analysis files.
        if system_changes is not None:
            self.remove_analysis()

        # Start writing the file.
        with open(filename, 'w') as f:
            # Write system name and label.
            f.write(format_fdf('SystemName', self.prefix))
            f.write(format_fdf('SystemLabel', self.prefix))
            f.write("\n")

            # Write explicitly given options first to
            # allow the user to override anything.
            fdf_arguments = self['fdf_arguments']
            keys = sorted(fdf_arguments.keys())
            for key in keys:
                f.write(format_fdf(key, fdf_arguments[key]))

            # Force siesta to return error on no convergence.
            # as default consistent with ASE expectations.
            if 'SCFMustConverge' not in fdf_arguments.keys():
                f.write(format_fdf('SCFMustConverge', True))
            f.write("\n")

            # Write spin level.
            f.write(format_fdf('Spin     ', self['spin']))
            # Spin backwards compatibility.
            if self['spin'] == 'collinear':
                f.write(format_fdf('SpinPolarized', (True, "# Backwards compatibility.")))
            elif self['spin'] == 'non-collinear':
                f.write(format_fdf('NonCollinear', (True, "# Backwards compatibility.")))



            # Write functional.
            functional, authors = self['xc']
            f.write(format_fdf('XC.functional', functional))
            f.write(format_fdf('XC.authors', authors))
            f.write("\n")

            # Write mesh cutoff and energy shift.
            f.write(format_fdf('MeshCutoff',
                               (self['mesh_cutoff'], 'eV')))
            f.write(format_fdf('PAO.EnergyShift',
                               (self['energy_shift'], 'eV')))
            f.write("\n")

            # Write the minimal arg
            self._write_species(f, atoms)
            self._write_structure(f, atoms)

            # Use the saved density matrix if only 'cell' and 'positions'
            # have changed.
            if (system_changes is None or
                ('numbers' not in system_changes and
                 'initial_magmoms' not in system_changes and
                 'initial_charges' not in system_changes)):
                f.write(format_fdf('DM.UseSaveDM', True))

            # Save density.
            if 'density' in properties:
                f.write(format_fdf('SaveRho', True))

            self._write_kpts(f)

            if self['bandpath'] is not None:
                lines = bandpath2bandpoints(self['bandpath'])
                f.write(lines)
                f.write('\n')

    def read(self, filename):
        """Read structural parameters from file .XV file
           Read other results from other files
           filename : siesta.XV
        """

        fname = self.getpath(filename)
        if not os.path.exists(fname):
            raise ReadError("The restart file '%s' does not exist" % fname)
        self.atoms = xv_to_atoms(fname)
        self.read_results()

    def getpath(self, fname=None, ext=None):
        """ Returns the directory/fname string """
        if fname is None:
            fname = self.prefix
        if ext is not None:
            fname = '{}.{}'.format(fname, ext)
        return os.path.join(self.directory, fname)

    def remove_analysis(self):
        """ Remove all analysis files"""
        filename = self.getpath(ext='RHO')
        if os.path.exists(filename):
            os.remove(filename)

    def _write_structure(self, f, atoms):
        """Translate the Atoms object to fdf-format.

        Parameters:
            - f:     An open file object.
            - atoms: An atoms object.
        """
        cell = atoms.cell
        f.write('\n')

        if cell.rank in [1, 2]:
            raise ValueError('Expected 3D unit cell or no unit cell.  You may '
                             'wish to add vacuum along some directions.')

        # Write lattice vectors
        if np.any(cell):
            f.write(format_fdf('LatticeConstant', '1.0 Ang'))
            f.write('%block LatticeVectors\n')
            for i in range(3):
                for j in range(3):
                    s = ('    %.15f' % cell[i, j]).rjust(16) + ' '
                    f.write(s)
                f.write('\n')
            f.write('%endblock LatticeVectors\n')
            f.write('\n')

        self._write_atomic_coordinates(f, atoms)

        # Write magnetic moments.
        magmoms = atoms.get_initial_magnetic_moments()

        # The DM.InitSpin block must be written to initialize to
        # no spin. SIESTA default is FM initialization, if the
        # block is not written, but  we must conform to the
        # atoms object.
        if self['spin'] != 'non-polarized':
            if sum(abs(magmoms)) == 0:
                f.write('#Empty block forces ASE initialization.\n')

            f.write('%block DM.InitSpin\n')
            for n, M in enumerate(magmoms):
                if M != 0:
                    f.write('    %d %.14f\n' % (n + 1, M))
            f.write('%endblock DM.InitSpin\n')
            f.write('\n')

    def _write_atomic_coordinates(self, f, atoms):
        """Write atomic coordinates.

        Parameters:
            - f:     An open file object.
            - atoms: An atoms object.
        """
        af = self.parameters.atomic_coord_format.lower()
        if af == 'xyz':
            self._write_atomic_coordinates_xyz(f, atoms)
        elif af == 'zmatrix':
            self._write_atomic_coordinates_zmatrix(f, atoms)
        else:
            raise RuntimeError('Unknown atomic_coord_format: {}'.format(af))

    def _write_atomic_coordinates_xyz(self, f, atoms):
        """Write atomic coordinates.

        Parameters:
            - f:     An open file object.
            - atoms: An atoms object.
        """
        species, species_numbers = self.species(atoms)
        f.write('\n')
        f.write('AtomicCoordinatesFormat  Ang\n')
        f.write('%block AtomicCoordinatesAndAtomicSpecies\n')
        for atom, number in zip(atoms, species_numbers):
            xyz = atom.position
            line = ('    %.9f' % xyz[0]).rjust(16) + ' '
            line += ('    %.9f' % xyz[1]).rjust(16) + ' '
            line += ('    %.9f' % xyz[2]).rjust(16) + ' '
            line += str(number) + '\n'
            f.write(line)
        f.write('%endblock AtomicCoordinatesAndAtomicSpecies\n')
        f.write('\n')

        origin = tuple(-atoms.get_celldisp().flatten())
        if any(origin):
            f.write('%block AtomicCoordinatesOrigin\n')
            f.write('     %.4f  %.4f  %.4f\n' % origin)
            f.write('%endblock AtomicCoordinatesOrigin\n')
            f.write('\n')

    def _write_atomic_coordinates_zmatrix(self, f, atoms):
        """Write atomic coordinates in Z-matrix format.

        Parameters:
            - f:     An open file object.
            - atoms: An atoms object.
        """
        species, species_numbers = self.species(atoms)
        f.write('\n')
        f.write('ZM.UnitsLength   Ang\n')
        f.write('%block Zmatrix\n')
        f.write('  cartesian\n')
        fstr = "{:5d}" + "{:20.10f}" * 3 + "{:3d}" * 3 + "{:7d} {:s}\n"
        a2constr = self.make_xyz_constraints(atoms)
        a2p, a2s = atoms.get_positions(), atoms.get_chemical_symbols()
        for ia, (sp, xyz, ccc, sym) in enumerate(zip(species_numbers,
                                                     a2p,
                                                     a2constr,
                                                     a2s)):
            f.write(fstr.format(
                sp, xyz[0], xyz[1], xyz[2], ccc[0],
                ccc[1], ccc[2], ia + 1, sym))
        f.write('%endblock Zmatrix\n')

        origin = tuple(-atoms.get_celldisp().flatten())
        if any(origin):
            f.write('%block AtomicCoordinatesOrigin\n')
            f.write('     %.4f  %.4f  %.4f\n' % origin)
            f.write('%endblock AtomicCoordinatesOrigin\n')
            f.write('\n')

    def make_xyz_constraints(self, atoms):
        """ Create coordinate-resolved list of constraints [natoms, 0:3]
        The elements of the list must be integers 0 or 1
          1 -- means that the coordinate will be updated during relaxation
          0 -- mains that the coordinate will be fixed during relaxation
        """
        from ase.constraints import FixAtoms, FixedLine, FixedPlane
        import sys
        import warnings

        a = atoms
        a2c = np.ones((len(a), 3), dtype=int)
        for c in a.constraints:
            if isinstance(c, FixAtoms):
                a2c[c.get_indices()] = 0
            elif isinstance(c, FixedLine):
                norm_dir = c.dir / np.linalg.norm(c.dir)
                if (max(norm_dir) - 1.0) > 1e-6:
                    raise RuntimeError(
                        'norm_dir: {} -- must be one of the Cartesian axes...'
                        .format(norm_dir))
                a2c[c.a] = norm_dir.round().astype(int)
            elif isinstance(c, FixedPlane):
                norm_dir = c.dir / np.linalg.norm(c.dir)
                if (max(norm_dir) - 1.0) > 1e-6:
                    raise RuntimeError(
                        'norm_dir: {} -- must be one of the Cartesian axes...'
                        .format(norm_dir))
                a2c[c.a] = abs(1 - norm_dir.round().astype(int))
            else:
                warnings.warn('Constraint {} is ignored at {}'
                              .format(str(c), sys._getframe().f_code))
        return a2c

    def _write_kpts(self, f):
        """Write kpts.

        Parameters:
            - f : Open filename.
        """
        if self["kpts"] is None:
            return
        kpts = np.array(self['kpts'])
        f.write('\n')
        f.write('#KPoint grid\n')
        f.write('%block kgrid_Monkhorst_Pack\n')

        for i in range(3):
            s = ''
            if i < len(kpts):
                number = kpts[i]
                displace = 0.0
            else:
                number = 1
                displace = 0
            for j in range(3):
                if j == i:
                    write_this = number
                else:
                    write_this = 0
                s += '     %d  ' % write_this
            s += '%1.1f\n' % displace
            f.write(s)
        f.write('%endblock kgrid_Monkhorst_Pack\n')
        f.write('\n')

    def _write_species(self, f, atoms):
        """Write input related the different species.

        Parameters:
            - f:     An open file object.
            - atoms: An atoms object.
        """
        species, species_numbers = self.species(atoms)

        if not self['pseudo_path'] is None:
            pseudo_path = self['pseudo_path']
        elif 'SIESTA_PP_PATH' in os.environ:
            pseudo_path = os.environ['SIESTA_PP_PATH']
        else:
            mess = "Please set the environment variable 'SIESTA_PP_PATH'"
            raise Exception(mess)

        f.write(format_fdf('NumberOfSpecies', len(species)))
        f.write(format_fdf('NumberOfAtoms', len(atoms)))

        pao_basis = []
        chemical_labels = []
        basis_sizes = []
        synth_blocks = []
        for species_number, spec in enumerate(species):
            species_number += 1
            symbol = spec['symbol']
            atomic_number = atomic_numbers[symbol]

            if spec['pseudopotential'] is None:
                if self.pseudo_qualifier() == '':
                    label = symbol
                    pseudopotential = label + '.psf'
                else:
                    label = '.'.join([symbol, self.pseudo_qualifier()])
                    pseudopotential = label + '.psf'
            else:
                pseudopotential = spec['pseudopotential']
                label = os.path.basename(pseudopotential)
                label = '.'.join(label.split('.')[:-1])

            if not os.path.isabs(pseudopotential):
                pseudopotential = join(pseudo_path, pseudopotential)

            if not os.path.exists(pseudopotential):
                mess = "Pseudopotential '%s' not found" % pseudopotential
                raise RuntimeError(mess)

            name = os.path.basename(pseudopotential)
            name = name.split('.')
            name.insert(-1, str(species_number))
            if spec['ghost']:
                name.insert(-1, 'ghost')
                atomic_number = -atomic_number

            name = '.'.join(name)
            pseudo_targetpath = self.getpath(name)

            if join(os.getcwd(), name) != pseudopotential:
                if islink(pseudo_targetpath) or isfile(pseudo_targetpath):
                    os.remove(pseudo_targetpath)
                symlink_pseudos = self['symlink_pseudos']

                if symlink_pseudos is None:
                    symlink_pseudos = not os.name == 'nt'

                if symlink_pseudos:
                    os.symlink(pseudopotential, pseudo_targetpath)
                else:
                    shutil.copy(pseudopotential, pseudo_targetpath)

            if not spec['excess_charge'] is None:
                atomic_number += 200
                n_atoms = sum(np.array(species_numbers) == species_number)

                paec = float(spec['excess_charge']) / n_atoms
                vc = get_valence_charge(pseudopotential)
                fraction = float(vc + paec) / vc
                pseudo_head = name[:-4]
                fractional_command = os.environ['SIESTA_UTIL_FRACTIONAL']
                cmd = '%s %s %.7f' % (fractional_command,
                                      pseudo_head,
                                      fraction)
                os.system(cmd)

                pseudo_head += '-Fraction-%.5f' % fraction
                synth_pseudo = pseudo_head + '.psf'
                synth_block_filename = pseudo_head + '.synth'
                os.remove(name)
                shutil.copyfile(synth_pseudo, name)
                synth_block = read_vca_synth_block(
                    synth_block_filename,
                    species_number=species_number)
                synth_blocks.append(synth_block)

            if len(synth_blocks) > 0:
                f.write(format_fdf('SyntheticAtoms', list(synth_blocks)))

            label = '.'.join(np.array(name.split('.'))[:-1])
            string = '    %d %d %s' % (species_number, atomic_number, label)
            chemical_labels.append(string)
            if isinstance(spec['basis_set'], PAOBasisBlock):
                pao_basis.append(spec['basis_set'].script(label))
            else:
                basis_sizes.append(("    " + label, spec['basis_set']))
        f.write((format_fdf('ChemicalSpecieslabel', chemical_labels)))
        f.write('\n')
        f.write((format_fdf('PAO.Basis', pao_basis)))
        f.write((format_fdf('PAO.BasisSizes', basis_sizes)))
        f.write('\n')

    def pseudo_qualifier(self):
        """Get the extra string used in the middle of the pseudopotential.
        The retrieved pseudopotential for a specific element will be
        'H.xxx.psf' for the element 'H' with qualifier 'xxx'. If qualifier
        is set to None then the qualifier is set to functional name.
        """
        if self['pseudo_qualifier'] is None:
            return self['xc'][0].lower()
        else:
            return self['pseudo_qualifier']

    def read_results(self):
        """Read the results.
        """
        self.read_number_of_grid_points()
        self.read_energy()
        self.read_forces_stress()
        self.read_eigenvalues()
        self.read_kpoints()
        self.read_dipole()
        self.read_pseudo_density()
        self.read_hsx()
        self.read_dim()
        if self.results['hsx'] is not None:
            self.read_pld(self.results['hsx'].norbitals,
                          len(self.atoms))
            self.atoms.cell = self.results['pld'].cell * Bohr
        else:
            self.results['pld'] = None

        self.read_wfsx()
        self.read_ion(self.atoms)

        self.read_bands()

    def read_bands(self):
        bandpath = self['bandpath']
        if bandpath is None:
            return

        if len(bandpath.kpts) < 1:
            return

        fname = self.getpath(ext='bands')
        with open(fname) as fd:
            kpts, energies, efermi = read_bands_file(fd)
        bs = resolve_band_structure(bandpath, kpts, energies, efermi)
        self.results['bandstructure'] = bs

    def band_structure(self):
        return self.results['bandstructure']

    def read_ion(self, atoms):
        """Read the ion.xml file of each specie
        """
        from ase.calculators.siesta.import_ion_xml import get_ion

        species, species_numbers = self.species(atoms)

        self.results['ion'] = {}
        for species_number, spec in enumerate(species):
            species_number += 1

            symbol = spec['symbol']
            atomic_number = atomic_numbers[symbol]

            if spec['pseudopotential'] is None:
                if self.pseudo_qualifier() == '':
                    label = symbol
                else:
                    label = '.'.join([symbol, self.pseudo_qualifier()])
                pseudopotential = self.getpath(label, 'psf')
            else:
                pseudopotential = spec['pseudopotential']
                label = os.path.basename(pseudopotential)
                label = '.'.join(label.split('.')[:-1])

            name = os.path.basename(pseudopotential)
            name = name.split('.')
            name.insert(-1, str(species_number))
            if spec['ghost']:
                name.insert(-1, 'ghost')
                atomic_number = -atomic_number
            name = '.'.join(name)

            label = '.'.join(np.array(name.split('.'))[:-1])

            if label not in self.results['ion']:
                fname = self.getpath(label, 'ion.xml')
                if os.path.isfile(fname):
                    self.results['ion'][label] = get_ion(fname)

    def read_hsx(self):
        """
        Read the siesta HSX file.
        return a namedtuple with the following arguments:
        'norbitals', 'norbitals_sc', 'nspin', 'nonzero',
        'is_gamma', 'sc_orb2uc_orb', 'row2nnzero', 'sparse_ind2column',
        'H_sparse', 'S_sparse', 'aB2RaB_sparse', 'total_elec_charge', 'temp'
        """
        from ase.calculators.siesta.import_functions import readHSX

        filename = self.getpath(ext='HSX')
        if isfile(filename):
            self.results['hsx'] = readHSX(filename)
        else:
            self.results['hsx'] = None

    def read_dim(self):
        """
        Read the siesta DIM file
        Retrun a namedtuple with the following arguments:
        'natoms_sc', 'norbitals_sc', 'norbitals', 'nspin',
        'nnonzero', 'natoms_interacting'
        """
        from ase.calculators.siesta.import_functions import readDIM

        filename = self.getpath(ext='DIM')
        if isfile(filename):
            self.results['dim'] = readDIM(filename)
        else:
            self.results['dim'] = None

    def read_pld(self, norb, natms):
        """
        Read the siesta PLD file
        Return a namedtuple with the following arguments:
        'max_rcut', 'orb2ao', 'orb2uorb', 'orb2occ', 'atm2sp',
        'atm2shift', 'coord_sc', 'cell', 'nunit_cells'
        """
        from ase.calculators.siesta.import_functions import readPLD

        filename = self.getpath(ext='PLD')
        if isfile(filename):
            self.results['pld'] = readPLD(filename, norb, natms)
        else:
            self.results['pld'] = None

    def read_wfsx(self):
        """
        Read the siesta WFSX file
        Return a namedtuple with the following arguments:
        """
        from ase.calculators.siesta.import_functions import readWFSX

        fname_woext = os.path.join(self.directory, self.prefix)

        if isfile(fname_woext + '.WFSX'):
            filename = fname_woext + '.WFSX'
            self.results['wfsx'] = readWFSX(filename)
        elif isfile(fname_woext + '.fullBZ.WFSX'):
            filename = fname_woext + '.fullBZ.WFSX'
            readWFSX(filename)
            self.results['wfsx'] = readWFSX(filename)
        else:
            self.results['wfsx'] = None

    def read_pseudo_density(self):
        """Read the density if it is there."""
        filename = self.getpath(ext='RHO')
        if isfile(filename):
            self.results['density'] = read_rho(filename)

    def read_number_of_grid_points(self):
        """Read number of grid points from SIESTA's text-output file. """

        fname = self.getpath(ext='out')
        with open(fname, 'r') as f:
            for line in f:
                line = line.strip().lower()
                if line.startswith('initmesh: mesh ='):
                    n_points = [int(word) for word in line.split()[3:8:2]]
                    self.results['n_grid_point'] = n_points
                    break
            else:
                raise RuntimeError

    def read_energy(self):
        """Read energy from SIESTA's text-output file.
        """
        fname = self.getpath(ext='out')
        with open(fname, 'r') as f:
            text = f.read().lower()

        assert 'final energy' in text
        lines = iter(text.split('\n'))

        # Get the energy and free energy the last time it appears
        for line in lines:
            has_energy = line.startswith('siesta: etot    =')
            if has_energy:
                self.results['energy'] = float(line.split()[-1])
                line = next(lines)
                self.results['free_energy'] = float(line.split()[-1])

        if ('energy' not in self.results or
            'free_energy' not in self.results):
            raise RuntimeError

    def read_forces_stress(self):
        """Read the forces and stress from the FORCE_STRESS file.
        """
        fname = self.getpath('FORCE_STRESS')
        with open(fname, 'r') as f:
            lines = f.readlines()

        stress_lines = lines[1:4]
        stress = np.empty((3, 3))
        for i in range(3):
            line = stress_lines[i].strip().split(' ')
            line = [s for s in line if len(s) > 0]
            stress[i] = [float(s) for s in line]

        self.results['stress'] = np.array(
            [stress[0, 0], stress[1, 1], stress[2, 2],
             stress[1, 2], stress[0, 2], stress[0, 1]])

        self.results['stress'] *= Ry / Bohr**3

        start = 5
        self.results['forces'] = np.zeros((len(lines) - start, 3), float)
        for i in range(start, len(lines)):
            line = [s for s in lines[i].strip().split(' ') if len(s) > 0]
            self.results['forces'][i - start] = [float(s) for s in line[2:5]]

        self.results['forces'] *= Ry / Bohr

    def read_eigenvalues(self):
        """ A robust procedure using the suggestion by Federico Marchesin """

        fname = self.getpath(ext='EIG')
        try:
            with open(fname, "r") as f:
                self.results['fermi_energy'] = float(f.readline())
                n, nspin, nkp = map(int, f.readline().split())
                _ee = np.split(
                    np.array(f.read().split()).astype(np.float), nkp)
        except (IOError):
            return 1

        ksn2e = np.delete(_ee, 0, 1).reshape([nkp, nspin, n])

        eigarray = np.empty((nspin, nkp, n))
        eigarray[:] = np.inf

        for k, sn2e in enumerate(ksn2e):
            for s, n2e in enumerate(sn2e):
                eigarray[s, k, :] = n2e

        assert np.isfinite(eigarray).all()

        self.results['eigenvalues'] = eigarray
        return 0

    def read_kpoints(self):
        """ Reader of the .KP files """

        fname = self.getpath(ext='KP')
        try:
            with open(fname, "r") as fd:
                nkp = int(next(fd))
                kpoints = np.empty((nkp, 3))
                kweights = np.empty(nkp)

                for i in range(nkp):
                    line = next(fd)
                    tokens = line.split()
                    numbers = np.array(tokens[1:]).astype(float)
                    kpoints[i] = numbers[:3]
                    kweights[i] = numbers[3]

        except (IOError):
            return 1

        self.results['kpoints'] = kpoints
        self.results['kweights'] = kweights

        return 0

    def read_dipole(self):
        """Read dipole moment. """
        dipole = np.zeros([1, 3])
        with open(self.getpath(ext='out'), 'r') as f:
            for line in f:
                if line.rfind('Electric dipole (Debye)') > -1:
                    dipole = np.array([float(f) for f in line.split()[5:8]])
        # debye to e*Ang
        self.results['dipole'] = dipole * 0.2081943482534

    def pyscf_tddft(
            self,
            Edir=np.array([1.0, 0.0, 0.0]),
            freq=np.arange(0.0, 10.0, 0.1),
            units='au',
            run_tddft=True,
            save_kernel=True,
            kernel_name="tddft_kernel.npy",
            fname="pol_tensor.npy",
            fname_nonin="noninpol_tensor.npy", **kw):
        """
        Perform TDDFT calculation using the pyscf.nao module for a molecule.

        Parameters
        ----------
        freq: array like
            frequency range for which the polarizability should
            be computed, in eV
        units : str, optional
            unit for the returned polarizability, can be au (atomic units)
            or nm**2
        run_tddft: to run the tddft_calculation or not
        fname: str
            Name of input file name for polariazbility tensor.
            if run_tddft is True: output file
            if run_tddft is False: input file

        kw: keywords for the tddft_iter function from pyscf

        Returns
        -------
            Add to the self.results dict the following items:
        freq range: array like
            array of dimension (nff) containing the frequency range in eV.

        polarizability nonin: array like (complex)
            array of dimension (nff, 3, 3) with nff the frequency number,
            the second and third dimension are the matrix elements of the
            non-interactive polarizability::

                P_xx, P_xy, P_xz, Pyx, .......


        polarizability: array like (complex)
            array of dimension (nff, 3, 3) with nff the frequency number,
            the second and third dimension are the matrix elements of the
            interactive polarizability::

                P_xx, P_xy, P_xz, Pyx, .......

        density change nonin: array like (complex)
            contains the non interacting density change in product basis

        density change inter: array like (complex)
            contains the interacting density change in product basis

        References
        ----------
        https://github.com/cfm-mpc/pyscf/tree/nao

        Example
        -------
        from ase.units import Ry, eV, Ha
        from ase.calculators.siesta import Siesta
        from ase import Atoms
        import numpy as np
        import matplotlib.pyplot as plt

        # Define the systems
        Na8 = Atoms('Na8',
                     positions=[[-1.90503810, 1.56107288, 0.00000000],
                                [1.90503810, 1.56107288, 0.00000000],
                                [1.90503810, -1.56107288, 0.00000000],
                                [-1.90503810, -1.56107288, 0.00000000],
                                [0.00000000, 0.00000000, 2.08495836],
                                [0.00000000, 0.00000000, -2.08495836],
                                [0.00000000, 3.22798122, 2.08495836],
                                [0.00000000, 3.22798122, -2.08495836]],
                     cell=[20, 20, 20])

        # Siesta input
        siesta = Siesta(
                    mesh_cutoff=150 * Ry,
                    basis_set='DZP',
                    pseudo_qualifier='',
                    energy_shift=(10 * 10**-3) * eV,
                    fdf_arguments={
                        'SCFMustConverge': False,
                        'COOP.Write': True,
                        'WriteDenchar': True,
                        'PAO.BasisType': 'split',
                        'DM.Tolerance': 1e-4,
                        'DM.MixingWeight': 0.01,
                        'MaxSCFIterations': 300,
                        'DM.NumberPulay': 4,
                        'XML.Write': True})

        Na8.set_calculator(siesta)
        e = Na8.get_potential_energy()
        freq, pol = siesta.get_polarizability_pyscf_inter(
            label="siesta",
            jcutoff=7,
            iter_broadening=0.15/Ha,
            xc_code='LDA,PZ',
            tol_loc=1e-6,
            tol_biloc=1e-7,
            freq = np.arange(0.0, 5.0, 0.05))
        # plot polarizability
        plt.plot(freq, pol[:, 0, 0].imag)
        plt.show()
        """

        from ase.calculators.siesta.mbpt_lcao_utils import pol2cross_sec
        assert units in ["nm**2", "au"]

        if run_tddft:
            from pyscf.nao import tddft_iter
            from ase.units import Ha

            tddft = tddft_iter(**kw)
            if save_kernel:
                np.save(kernel_name, tddft.kernel)

            omegas = freq / Ha + 1j * tddft.eps
            tddft.comp_dens_nonin_along_Eext(omegas, Eext=Edir)
            tddft.comp_dens_inter_along_Eext(omegas, Eext=Edir)

            # save polarizability tensor and density change to files
            self.results["freq range"] = freq
            self.results['polarizability nonin'] = np.zeros(
                (freq.size, 3, 3),
                dtype=tddft.p0_mat.dtype)
            self.results['polarizability inter'] = np.zeros(
                (freq.size, 3, 3),
                dtype=tddft.p_mat.dtype)
            self.results["density change nonin"] = tddft.dn0
            self.results["density change inter"] = tddft.dn
            for xyz1 in range(3):
                for xyz2 in range(3):
                    if units == 'nm**2':
                        p0 = pol2cross_sec(-tddft.p0_mat[xyz1, xyz2, :], freq)
                        p = pol2cross_sec(-tddft.p_mat[xyz1, xyz2, :], freq)
                        self.results['polarizability nonin'][:, xyz1, xyz2] = \
                            p0
                        self.results['polarizability inter'][:, xyz1, xyz2] = \
                            p
                    else:
                        self.results['polarizability nonin'][:, xyz1, xyz2] = \
                            -tddft.p0_mat[xyz1, xyz2, :]
                        self.results['polarizability inter'][:, xyz1, xyz2] = \
                            -tddft.p_mat[xyz1, xyz2, :]

        else:
            # load polarizability tensor from previous calculations
            p0_mat = np.load(fname_nonin)
            p_mat = np.load(fname)

            self.results['polarizability nonin'] = np.zeros(
                (freq.size, 3, 3),
                dtype=p0_mat.dtype)
            self.results['polarizability inter'] = np.zeros(
                (freq.size, 3, 3),
                dtype=p_mat.dtype)

            for xyz1 in range(3):
                for xyz2 in range(3):
                    if units == 'nm**2':
                        p0 = pol2cross_sec(-p0_mat[xyz1, xyz2, :], freq)
                        p = pol2cross_sec(-p_mat[xyz1, xyz2, :], freq)

                        self.results['polarizability nonin'][:, xyz1, xyz2] = \
                            p0
                        self.results['polarizability inter'][:, xyz1, xyz2] = \
                            p
                    else:
                        self.results['polarizability nonin'][:, xyz1, xyz2] = \
                            -p0_mat[xyz1, xyz2, :]
                        self.results['polarizability inter'][:, xyz1, xyz2] = \
                            -p_mat[xyz1, xyz2, :]

    def pyscf_tddft_eels(
            self,
            velec=np.array([20.0, 0.0, 0.0]),
            b=np.array([0.0, 0.0, 0.0]),
            freq=np.arange(0.0, 10.0, 0.1),
            tddft=None,
            save_kernel=True,
            kernel_name="tddft_kernel.npy",
            tmp_fname=None,
            **kw):
        r"""
        Perform TDDFT calculation using the pyscf.nao module for a molecule.
        The external pertubation is created by a electron moving at
        the velocity velec and with an impact parameter b.

        Parameters
        ----------
        freq: array like
            frequency range for which the polarizability should
            be computed, in eV
        velec: array like
            velocity vector of the projectile
        b: array like
            offset vector of the projectile
        tddft: tddft_tem class from a previous calculation
        save_kernel: save the kernel for future use
        kernel_name: name of the file for the kernel
        tmp_fname: temporary name to save the eels spectra
        while running the calculations.
        kw: keywords for the tddft_tem function from pyscf

        Returns
        -------
        tddft:
            if running pyscf_tddft_eels in a loop over the velocity or the
            impact parameter, there is no point to initialize again the tddft
            calculation (vertex and kernel will be the same)

            Add to the self.results dict the following items:
        freq range: array like
            array of dimension (nff) containing the frequency range in eV.

        eel spectra nonin: array like (complex)
            array of dimension (nff) with nff the frequency number,


        eel spectra inter: array like (complex)
            array of dimension (nff) with nff the frequency number,

        density change eels nonin: array like (complex)
            contains the non interacting density change in product basis

        density change eels inter: array like (complex)
            contains the interacting density change in product basis

        References
        ----------
        https://github.com/cfm-mpc/pyscf/tree/nao

        Example
        -------
        from ase.units import Ry, eV, Ha
        from ase.calculators.siesta import Siesta
        from ase import Atoms
        import numpy as np
        import matplotlib.pyplot as plt

        # Define the systems
        Na8 = Atoms('Na8',
                    positions=[[-1.90503810, 1.56107288, 0.00000000],
                               [1.90503810, 1.56107288, 0.00000000],
                               [1.90503810, -1.56107288, 0.00000000],
                               [-1.90503810, -1.56107288, 0.00000000],
                               [0.00000000, 0.00000000, 2.08495836],
                               [0.00000000, 0.00000000, -2.08495836],
                               [0.00000000, 3.22798122, 2.08495836],
                               [0.00000000, 3.22798122, -2.08495836]],
                    cell=[20, 20, 20])

        # enter siesta input
        siesta = Siesta(
            mesh_cutoff=150 * Ry,
            basis_set='DZP',
            pseudo_qualifier='',
            energy_shift=(10 * 10**-3) * eV,
            fdf_arguments={
                'SCFMustConverge': False,
                'COOP.Write': True,
                'WriteDenchar': True,
                'PAO.BasisType': 'split',
                'DM.Tolerance': 1e-4,
                'DM.MixingWeight': 0.01,
                'MaxSCFIterations': 300,
                'DM.NumberPulay': 4,
                'XML.Write': True})


        Na8.set_calculator(siesta)
        e = Na8.get_potential_energy()
        tddft = siesta.pyscf_tddft_eels(
            label="siesta", jcutoff=7, iter_broadening=0.15/Ha,
            xc_code='LDA,PZ', tol_loc=1e-6, tol_biloc=1e-7,
            freq=np.arange(0.0, 5.0, 0.05))

        # plot eel spectra
        fig = plt.figure(1)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.plot(siesta.results["freq range"],
                 siesta.results["eel spectra nonin"].imag)
        ax2.plot(siesta.results["freq range"],
                 siesta.results["eel spectra inter"].imag)

        ax1.set_xlabel(r"$\omega$ (eV)")
        ax2.set_xlabel(r"$\omega$ (eV)")

        ax1.set_ylabel(r"Im($P_{xx}$) (au)")
        ax2.set_ylabel(r"Im($P_{xx}$) (au)")

        ax1.set_title(r"Non interacting")
        ax2.set_title(r"Interacting")

        fig.tight_layout()
        plt.show()
        """
        from pyscf.nao import tddft_tem
        from ase.units import Ha

        assert velec.size == 3
        assert b.size == 3

        if tddft is None:
            self.results["freq range"] = freq
            omegas = freq / Ha

            # for eels, omega is real array
            tddft = tddft_tem(freq=omegas, **kw)
            if save_kernel:
                np.save(kernel_name, tddft.kernel)

        self.results['eel spectra nonin'] = \
            tddft.get_spectrum_nonin(velec=velec,
                                     beam_offset=b,
                                     tmp_fname=tmp_fname)

        self.results['eel spectra inter'] = \
            tddft.get_spectrum_inter(velec=velec,
                                     beam_offset=b,
                                     tmp_fname=tmp_fname)

        self.results["density change eels nonin"] = tddft.dn0
        self.results["density change eels inter"] = tddft.dn

        return tddft

    def get_polarizability_mbpt(self, mbpt_inp=None,
                                output_name='mbpt_lcao.out',
                                format_output='hdf5', units='au'):
        """
        Warning!!
            Out dated version, try get_polarizability_pyscf


        Calculate the polarizability by running the mbpt_lcao program.
        The mbpt_lcao program need the siesta output, therefore siesta need
        to be run first.

        Parameters
        ----------
        mbpt_inp : dict, optional
            dictionnary of the input for the mbpt_lcao program
            (http://mbpt-domiprod.wikidot.com/list-of-parameters)
            if mbpt_inp is None, the function read the output file
            from a previous mbpt_lcao run.
        output_name : str, optional
            Name of the mbpt_lcao output
        format_output : str, optional
            Format of the mbpt_lcao output data,
            if hdf5, the output name is tddft_iter_output.hdf5 if
            do_tddft_iter is set to 1 the output name is
            tddft_tem_output.hdf5 if do_tddft_tem is set to 1
            if txt, a lot of output data files are produced depending on
            the input, in the text and fortran binaries format
        units : str, optional
            unit for the returned polarizability, can be au (atomic units)
            or nm**2

        Returns
        -------
        freq : array like
            array of dimension (nff) containing the frequency range in eV.

        self.results['polarizability'], array like
            array of dimension (nff, 3, 3) with nff the frequency number,
            the second and third dimension are the matrix elements of the
            polarizability::

                P_xx, P_xy, P_xz, Pyx, .......

        References
        ----------
        http://mbpt-domiprod.wikidot.com

        Example
        -------
        import os
        from ase.units import Ry, eV
        from ase.calculators.siesta import Siesta
        from ase import Atoms
        import numpy as np
        import matplotlib.pyplot as plt

        #Define the systems
        Na8 = Atoms('Na8',
        positions=[[-1.90503810, 1.56107288, 0.00000000],
                    [1.90503810, 1.56107288, 0.00000000],
                    [1.90503810, -1.56107288, 0.00000000],
                    [-1.90503810, -1.56107288, 0.00000000],
                    [0.00000000, 0.00000000, 2.08495836],
                    [0.00000000, 0.00000000, -2.08495836],
                    [0.00000000, 3.22798122, 2.08495836],
                    [0.00000000, 3.22798122, -2.08495836]],
                    cell=[20, 20, 20])

        #enter siesta input
        siesta = Siesta(
            mesh_cutoff=150 * Ry,
            basis_set='DZP',
            pseudo_qualifier='',
            energy_shift=(10 * 10**-3) * eV,
            fdf_arguments={
            'SCFMustConverge': False,
            'COOP.Write': True,
            'WriteDenchar': True,
            'PAO.BasisType': 'split',
            'DM.Tolerance': 1e-4,
            'DM.MixingWeight': 0.01,
            'MaxSCFIterations': 300,
            'DM.NumberPulay': 4})


        #mbpt_lcao input
        mbpt_inp = {'prod_basis_type' : 'MIXED',
                    'solver_type' : 1,
                    'gmres_eps' : 0.001,
                    'gmres_itermax':256,
                    'gmres_restart':250,
                    'gmres_verbose':20,
                    'xc_ord_lebedev':14,
                    'xc_ord_gl':48,
                    'nr':512,
                    'akmx':100,
                    'eigmin_local':1e-06,
                    'eigmin_bilocal':1e-08,
                    'freq_eps_win1':0.15,
                    'd_omega_win1':0.05,
                    'dt':0.1,
                    'omega_max_win1':5.0,
                    'ext_field_direction':2,
                    'dr':np.array([0.3, 0.3, 0.3]),
                    'para_type':'MATRIX',
                    'chi0_v_algorithm':14,
                    'format_output':'text',
                    'comp_dens_chng_and_polarizability':1,
                    'store_dens_chng':1,
                    'enh_given_volume_and_freq':0,
                    'diag_hs':0,
                    'do_tddft_tem':0,
                    'do_tddft_iter':1,
                    'plot_freq':3.02,
                    'gwa_initialization':'SIESTA_PB'}


        Na8.set_calculator(siesta)
        e = Na8.get_potential_energy() #run siesta
        freq, pol = siesta.get_polarizability_siesta(mbpt_inp,
                                                     format_output='txt',
                                                     units='nm**2')

        #plot polarizability
        plt.plot(freq, pol[:, 0, 0])

        plt.show()
        """
        from ase.calculators.siesta.mbpt_lcao import MBPT_LCAO
        from ase.calculators.siesta.mbpt_lcao_io import read_mbpt_lcao_output

        warnings.warn("Out dated version, try get_polarizability_pyscf")

        if mbpt_inp is not None:
            tddft = MBPT_LCAO(mbpt_inp)
            tddft.run_mbpt_lcao(output_name, True)

        r = read_mbpt_lcao_output()

        r.args.format_input = format_output

        # read real part
        r.args.ReIm = 're'
        data = r.Read()
        self.results['polarizability'] = data.Array

        # read imaginary part
        r.args.ReIm = 'im'
        data = r.Read()
        self.results['polarizability'] = (self.results['polarizability'] +
                                          complex(0.0, 1.0) * data.Array)

        if units == 'nm**2':
            from ase.calculators.siesta.mbpt_lcao_utils import pol2cross_sec
            for i in range(2):
                for j in range(2):
                    p = pol2cross_sec(self.results['polarizability'][:, i, j],
                                      data.freq)
                    self.results['polarizability'][:, i, j] = p

            print('unit nm**2')
            # self.results['polarizability'] = data.Array
        elif units == 'au':
            print('unit au')
            # self.results['polarizability'] = data.Array
        else:
            raise ValueError('units can be only au or nm**2')

        return data.freq, self.results['polarizability']

    def get_fermi_level(self):
        return self.results['fermi_energy']

    def get_k_point_weights(self):
        return self.results['kweights']

    def get_ibz_k_points(self):
        return self.results['kpoints']


class Siesta3_2(Siesta):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            "The Siesta3_2 calculator class will no longer be supported. "
            "Use 'ase.calculators.siesta.Siesta in stead. "
            "If using the ASE interface with SIESTA 3.2 you must explicitly "
            "include the keywords 'SpinPolarized', 'NonCollinearSpin' and "
            "'SpinOrbit' if needed.",
            np.VisibleDeprecationWarning)
        Siesta.__init__(self, *args, **kwargs)
