# -*- coding: utf-8 -*-
"""This module defines an interface to ONETEP for use by the ASE.

Authors:
    Edward Tait, ewt23@cam.ac.uk
    Nicholas Hine, n.d.m.hine@warwick.ac.uk (current maintainer)

    Based on castep.py by:
    Max Hoffmann, max.hoffmann@ch.tum.de
    Jörg Meyer, joerg.meyer@ch.tum.de
"""

from copy import deepcopy
from os.path import isfile
from warnings import warn

from numpy import array

from ase import Atoms
from ase.calculators.calculator import FileIOCalculator, ReadError
from ase.parallel import paropen
from ase.units import Bohr, Hartree


__all__ = ['Onetep']


class Onetep(FileIOCalculator):
    """Implements the calculator for the onetep linear
    scaling DFT code. Recommended ASE_ONETEP_COMMAND format
    is "onetep_executable_name PREFIX.dat > PREFIX.out 2> PREFIX.err" """

    implemented_properties = ['energy', 'forces', 'dipole', 'magmom']

    # Used to indicate 'parameters' which shouldn't be written to
    # the onetep input file in the standard <key> : <value> format
    # for example the NGWF radius is used in the species block and isn't
    # written elsewhere in the input file
    _dummy_parameters = ['ngwf_radius', 'xc', 'species_ngwf_radius',
                         'species_ngwf_number', 'species_solver',
                         'ngwf_radius_cond', 'pseudo_suffix',
                         'species_pseudo', 'species_core_wf',
                         'species_solver_cond', 'species_ngwf_number_cond',
                         'species_ngwf_radius_cond']

    # Used to indicate which parameters are a kpoint path and should be
    # written as such
    _path_parameters = ['bsunfld_kpoint_path', 'bs_kpoint_path']

    # Used to indicate which parameters are a block listing atom
    # groupings for a variety of purposes
    _group_parameters = ['species_bsunfld_groups', 'species_ldos_groups',
                         'species_locdipole_groups',
                         'species_bsunfld_projatoms',
                         'species_pdos_groups', 'species_tddft_ct',
                         'species_tddft_kernel', 'nbo_write_species',
                         'species_ngwf_plot']

    # Used to indicate which parameters are a block of any other sort
    # other than those above (the contents of the parameter is reproduced
    # verbatim within the block)
    _block_parameters = _path_parameters + _group_parameters + [
                        'species_constraints', 'nbo_species_ngwflabel',
                        'ddec_rmse_vdw', 'vdw_params', 'sol_ions', 'swri']

    default_parameters = {'cutoff_energy': '1000 eV',
                          'kernel_cutoff': '1000 bohr',
                          'ngwf_radius': 12.0,
                          'ngwf_radius_cond': -1.0}

    name = 'onetep'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, command=None, atoms=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command, **kwargs)

        self.species = []
        self.species_cond = []
        self.pseudos = []
        self.core_wfs = []
        self.solvers = []
        self.solvers_cond = []
        self.restart = False
        self.prefix = label
        self.directory = '.'

    def read(self, label):
        """Read a onetep .out file into the current instance."""

        FileIOCalculator.read(self, label)

        onetep_file = self.label + '.out'

        warnings = []

        try:
            out = paropen(onetep_file, 'r')
        except IOError:
            raise ReadError('Could not open output file "%s"' % onetep_file)

        # keep track of what we've read in
        read_lattice = False
        read_species = False
        read_positions = False

        line = out.readline()

        if self.atoms is None:
            self.atoms = Atoms()
            self.atoms.calc = self

        while line:
            clean_line = line.strip().lower()
            if '%block lattice_cart' in clean_line:
                self._read_lattice(out)
                read_lattice = True
            elif '%block species_pot' in clean_line:
                self._read_species_pot(out)
            elif clean_line.endswith('%block species_atomic_set'):
                self._read_species_solver(out)
            elif clean_line.endswith('%block species'):
                self._read_species(out)
                read_species = True
            elif '%block positions_abs' in clean_line:
                self._read_positions(out)
                read_positions = True
            elif '%block species_cond' in clean_line:
                self._read_species(out,cond=True)
            elif '%block species_atomic_set_cond' in clean_line:
                self._read_species_solver(out,cond=True)
            elif 'warn' in line.lower():
                warnings.append(line)
            line = out.readline()
        out.close()

        if warnings:
            warn('WARNING: %s contains warnings' % onetep_file)
            for warning in warnings:
                warn(warning)

        if not (read_lattice and read_species and read_positions):
            raise ReadError('Failed to read in essential calculation'
                            ' data from output file "%s"' % onetep_file)

        self.read_results(label)

    def read_results(self,label=None):
        FileIOCalculator.read_results(self)

        if label is None:
           onetep_file = self.label + '.out'
        else:
           onetep_file = label + '.out'

        warnings = []

        try:
            out = paropen(onetep_file, 'r')
        except IOError:
            raise ReadError('Could not open output file "%s"' % onetep_file)

        line = out.readline()
        while line:
            if '| Total' in line:
                self.results['energy'] = Hartree * float(line.split()[-2])
            elif ('Element  Atom         Cartesian components (Eh/a)'
                  in line):
                self._read_forces(out)
            elif ('Final Configuration' in line):
                self._read_geom_output(out)
            elif ('Integrated spin density' in line):
                self.results['magmom'] = self._read_magmom(line)
            elif '|Excitation|    Energy (in Ha)   |     Oscillator Str' in line:
                self._read_excitations(out)
            elif ('Dipole Moment Calculation' in line):
                self.results['dipole'] = self._read_dipole(out)
            elif 'warn' in line.lower():
                warnings.append(line)
            line = out.readline()

        if warnings:
            warn('WARNING: %s contains warnings' % onetep_file)
            for warning in warnings:
                warn(warning)

    def _read_lattice(self, out):
        """ read the lattice parameters out of a onetep .out formatted file
        stream"""

        axes = []

        l = out.readline()
        # onetep assumes lengths are in atomic units by default
        conv_fac = Bohr
        if 'ang' in l:
            l = out.readline()
            conv_fac = 1.0
        elif 'bohr' in l:
            l = out.readline()

        for _ in range(0, 3):
            l = l.strip()
            p = l.split()
            if len(p) != 3:
                raise ReadError('Malformed Lattice block line "%s"' % l)
            try:
                axes.append([conv_fac * float(comp) for comp in p[0:3]])
            except ValueError:
                raise ReadError("Can't parse line \"%s\" in axes block" % l)
            l = out.readline()
        self.atoms.set_cell(axes)

    def _read_positions(self, out):
        """Read the contents of a positions_abs block into the calculator's
        atoms object, setting both species and positions. Tries to strip out
        comment lines and is aware of angstom vs. bohr"""

        line = out.readline()
        # onetep assumes lengths are in atomic units by default
        conv_fac = Bohr
        if 'ang' in line:
            line = out.readline()
            conv_fac = 1.0
        elif 'bohr' in line:
            line = out.readline()
        symbols = []
        positions = []
        while '%endblock' not in line.lower():
            line = line.strip()
            if line[0] != '#':
                atom, suffix = line.split(None, 1)
                pos = suffix.split(None, 3)[0:3]
                try:
                    pos = [conv_fac * float(p) for p in pos]
                except ValueError:
                    raise ReadError('Malformed position line "%s"', line)
                symbols.append(atom)
                positions.append(pos)
            line = out.readline()
        tags = deepcopy(symbols)
        for j in range(len(symbols)):
            symbols[j] = ''.join(i for i in symbols[j] if not i.isdigit())
        for j in range(len(tags)):
            tags[j] = ''.join(i for i in tags[j] if not i.isalpha())
            if tags[j]=='':
               tags[j]='0'
            tags[j] = int(tags[j])
        if len(self.atoms)!=len(symbols):
           self.atoms = Atoms(symbols=symbols,positions=positions)
        self.atoms.set_chemical_symbols(symbols)
        self.atoms.set_tags(tags)
        self.atoms.set_positions(positions)

    def _read_dipole(self, out):
        """Reads total dipole moment from ONETEP output file"""

        # Find start of total dipole moment block
        line = ()
        while 'Total dipole moment' not in line:
            line = out.readline()

        # Read total dipole moment
        dipolemoment = []
        for label, pos in sorted({'dx': 6, 'dy': 2, 'dz': 2}.items()):
            assert label in line.split()
            value = float(line.split()[pos])*Bohr
            dipolemoment.append(value)
            line = out.readline()

        return array(dipolemoment)

    def _read_magmom(self, line):
        """Reads magnetic moment from Integrated Spin line"""
        return float(line.split()[4])

    def _read_geom_output(self, out):
        """Reads geometry optimisation output from ONETEP output file"""
        conv_fac = Bohr

        # Find start of atom positions
        while 'x-----' not in out.readline():
            pass
        symbols = []
        positions = []
        # Read atom positions
        line = out.readline()
        while 'xxxxxx' not in line:
            line = line.strip()
            pos = line.split()[3:6]
            pos = [conv_fac * float(p) for p in pos]
            atom = line.split()[1]
            positions.append(pos)
            symbols.append(atom)
            line = out.readline()
        if len(positions) != len(self.atoms):
            raise ReadError('Wrong number of atoms found in output geometry'
                            'block')
        if len(symbols) != len(self.atoms):
            raise ReadError('Wrong number of atoms found in output geometry'
                            'block')

        # Update atoms object with new positions (and symbols)
        self.atoms.set_positions(positions)
        self.atoms.set_chemical_symbols(symbols)

    def _read_species(self, out, cond=False):
        """ Read in species block from a onetep output file"""
        line = out.readline().strip()
        species = []
        while '%endblock' not in line.lower():
            atom, element, z, nngwf, ngwf_radius = line.split(None, 5)
            z = int(z)
            nngwf = int(nngwf)
            ngwf_radius = float(ngwf_radius)
            species.append((atom, element, z, nngwf, ngwf_radius,))
            line = out.readline().strip()
        if not cond:
            self.set_species(species)
        else:
            self.set_species_cond(species)

    def _read_species_pot(self, out):
        """ Read in pseudopotential information from a onetep output file"""
        line = out.readline().strip()
        pots = []
        while '%endblock' not in line.lower() and len(line) > 0:
            atom, suffix = line.split(None, 1)
            filename = suffix.split('#', 1)[0].strip()
            filename = filename.replace('"', '')   # take out quotes
            filename = filename.replace("'", '')
            pots.append((atom, filename,))
            line = out.readline().strip()
        if len(line) == 0:
            raise ReadError('End of file while reading potential block')
        self.set_pseudos(pots)

    def _read_species_solver(self, out, cond=False):
        """ Read in pseudopotential information from a onetep output file"""
        line = out.readline().strip()
        solvers = []
        while '%endblock' not in line.lower() and len(line) > 0:
            atom, suffix = line.split(None, 1)
            solver_str = suffix.split('#', 1)[0].strip()
            solvers.append((atom, solver_str))
            line = out.readline().strip()
        if len(line) == 0:
            raise ReadError('End of file while reading solver block')
        if not cond:
            self.set_solvers(solvers)
        else:
            self.set_solvers_cond(solvers)

    def _read_forces(self, out):
        """ Extract the computed forces from a onetep output file"""
        forces = []
        atomic2ang = Hartree / Bohr
        while True:
            line = out.readline()
            fields = line.split()
            if len(fields) > 6:
                break
        while len(fields) == 7:
            force = [float(fcomp) * atomic2ang for fcomp in fields[-4:-1]]
            forces.append(force)
            line = out.readline()
            fields = line.split()
        self.results['forces'] = array(forces)

    def _read_excitations(self,out):
        """ Extract the computed electronic excitations from a onetep output
        file."""
        excitations = []
        line = out.readline()
        while line:
            words = line.split()
            if len(words)==0:
                break
            excitations.append([float(words[0]),float(words[1])*Hartree,float(words[2])])
            line = out.readline()
        self.results['excitations'] = array(excitations)

    def _generate_species_block(self, cond=False):
        """Create a default onetep species block, use -1 for the NGWF number
        to trigger automatic NGWF number assigment using onetep's internal
        routines."""

        # check if we need to do anything.
        if len(self.species) == len(self.atoms.get_chemical_symbols()):
            return

        parameters = self.parameters

        atoms = self.atoms
        if not cond:
            self.species = []
            default_ngwf_radius = self.parameters['ngwf_radius']
            species_ngwf_rad_var = 'species_ngwf_radius'
            species_ngwf_num_var = 'species_ngwf_number'
        else:
            self.species_cond = []
            default_ngwf_radius = self.parameters['ngwf_radius_cond']
            species_ngwf_rad_var = 'species_ngwf_radius_cond'
            species_ngwf_num_var = 'species_ngwf_number_cond'
        for sp in set(zip(atoms.get_atomic_numbers(),
                          atoms.get_chemical_symbols(),
                          ["" if i==0 else str(i) for i in atoms.get_tags()])):
            try:
                ngrad = parameters[species_ngwf_rad_var][sp[1]]
            except KeyError:
                ngrad = default_ngwf_radius
            try:
                ngnum = parameters[species_ngwf_num_var][sp[1]]
            except KeyError:
                ngnum = -1
            if not cond:
                self.species.append((sp[1]+sp[2], sp[1], sp[0], ngnum, ngrad))
            else:
                self.species_cond.append((sp[1]+sp[2], sp[1], sp[0], ngnum, ngrad))

    def _generate_pseudo_block(self):
        """Create a default onetep pseudopotentials block, using the
        element name with the variable pseudo_suffix appended to it by default,
        unless the user has set overrides for specific species by setting
        specific entries in species_pseudo"""

        for sp in self.species:
            try:
                pseudo_string = self.parameters['species_pseudo'][sp[0]]
            except KeyError:
                try:
                    pseudo_string = sp[1] + self.parameters['pseudo_suffix']
                except KeyError:
                    pseudo_string = sp[1] # bare elem name if pseudo suffix empty
            self.pseudos.append((sp[0], pseudo_string))

    def _generate_solver_block(self,cond=False):
        """Create a default onetep pseudoatomic solvers block, using 'SOLVE'
        unless the user has set overrides for specific species by setting
        specific entries in species_solver (_cond)"""

        if not cond:
            solver_var = 'species_solver'
        else:
            solver_var = 'species_solver_cond'
        for sp in self.species:
            try:
                atomic_string = self.parameters[solver_var][sp[0]]
            except KeyError:
                atomic_string = 'SOLVE'
            if not cond:
                self.solvers.append((sp[0], atomic_string))
            else:
                self.solvers_cond.append((sp[0],atomic_string))

    def _generate_core_wf_block(self):
        """Create a default onetep core wavefunctions block, using 'NONE'
        unless the user has set overrides for specific species by setting
        specific entries in species_core_wf. If all are NONE, no block
        will be printed"""

        any_core_wfs = False
        for sp in self.species:
            try:
                core_wf_string = self.parameters['species_core_wf'][sp[0]]
                any_core_wfs = True
            except KeyError:
                core_wf_string = 'NONE'

            self.core_wfs.append((sp[0], core_wf_string))

        # if no species core wavefunction definitions were set to anything
        # other than 'NONE', delete the block entirely
        if not any_core_wfs:
            self.core_wfs = []

    def set_pseudos(self, pots):
        """ Sets the pseudopotential files used in this dat file """
        self.pseudos = deepcopy(pots)

    def set_solvers(self, solvers):
        """ Sets the solver strings used in this dat file """
        self.solvers = deepcopy(solvers)

    def set_solvers_cond(self, solvers):
        """ Sets the solver strings used in this dat file """
        self.solvers_cond = deepcopy(solvers)

    def set_atoms(self, atoms):
        self.atoms = atoms

    def set_species(self, sp):
        """ Sets the species in the current dat instance,
        in onetep this includes both atomic number information
        as well as NGWF parameters like number and cut off radius"""
        self.species = deepcopy(sp)

    def set_species_cond(self, spc):
        """ Sets the conduction species in the current dat instance,
        in onetep this includes both atomic number information
        as well as NGWF parameters like number and cut off radius"""
        self.species_cond = deepcopy(spc)

    def write_input(self, atoms, properties=None, system_changes=None):
        """Only writes the input .dat file and return
        This can be useful if one quickly needs to prepare input files
        for a cluster where no python or ASE is available. One can than
        upload the file manually and read out the results using
        Onetep().read().
        """

        if atoms is None:
            atoms = self.atoms

        if self.restart:
            self.parameters['read_tightbox_ngwfs'] = True
            self.parameters['read_denskern'] = True

        self._generate_species_block()

        if len(self.pseudos) < len(self.species):
            if 'pseudo_suffix' in self.parameters:
                self._generate_pseudo_block()

        if len(self.solvers) < len(self.species):
            self._generate_solver_block()

        if 'ngwf_radius_cond' in self.parameters:
            if len(self.species_cond) < len(self.species):
                self._generate_species_block(cond=True)
            if len(self.solvers_cond) < len(self.species):
                self._generate_solver_block(cond=True)

        if len(self.core_wfs) < len(self.species):
            self._generate_core_wf_block()

        self._write_dat()

    def get_dipole_moment(self, atoms=None):
        self.parameters['polarisation_calculate'] = True
        self.parameters['do_properties'] = True
        return FileIOCalculator.get_dipole_moment(self, atoms)

    def get_forces(self, atoms=None):
        self.parameters['write_forces'] = True
        return FileIOCalculator.get_forces(self, atoms)

    def _write_dat(self, force_write=True):
        """This export function write minimal information to
        a .dat file. If the atoms object is a trajectory, it will
        take the last image.
        """
        filename = self.label + '.dat'

        if self.atoms is None:
            raise Exception('No associated atoms object.')

        atoms = self.atoms
        parameters = self.parameters

        if isfile(filename) and not force_write:
            raise Exception('Target input file already exists.')

        if 'xc' in parameters and 'xc_functional' in parameters \
                and parameters['xc'] != parameters['xc_functional']:
            raise Exception('Conflicting functionals defined! %s vs. %s' %
                            (parameters['xc'], parameters['xc_functional']))

        fd = open(filename, 'w')
        fd.write('######################################################\n')
        fd.write('#ONETEP .dat file: %s\n' % filename)
        fd.write('#Created using the Atomic Simulation Environment (ASE)\n')
        fd.write('######################################################\n\n')

        fd.write('%BLOCK LATTICE_CART\n')
        fd.write('ang\n')
        for line in atoms.get_cell():
            fd.write('    %.10f %.10f %.10f\n' % tuple(line))
        fd.write('%ENDBLOCK LATTICE_CART\n\n\n')

        keyword = 'POSITIONS_ABS'
        positions = atoms.get_positions()
        tags = ["" if i==0 else str(i) for i in atoms.get_tags()]
        pos_block = [('%s %8.6f %8.6f %8.6f' %
                      (x+z, y[0], y[1], y[2])) for (x, y, z)
                     in zip(atoms.get_chemical_symbols(), positions, tags)]

        fd.write('%%BLOCK %s\n' % keyword)
        fd.write('ang\n')
        for line in pos_block:
            fd.write('    %s\n' % line)
        fd.write('%%ENDBLOCK %s\n\n' % keyword)

        keyword = 'SPECIES'
        sp_block = [('%s %s %d %d %8.6f' % sp) for sp in self.species]
        fd.write('%%BLOCK %s\n' % keyword)
        for line in sorted(sp_block):
            fd.write('    %s\n' % line)
        fd.write('%%ENDBLOCK %s\n\n' % keyword)

        if ((self.parameters['ngwf_radius_cond'] > 0) or
            len(self.species_cond)==len(self.species)):
            keyword = 'SPECIES_COND'
            sp_block = [('%s %s %d %d %8.6f' % sp) for sp in self.species_cond]
            fd.write('%%BLOCK %s\n' % keyword)
            for line in sorted(sp_block):
                fd.write('    %s\n' % line)
            fd.write('%%ENDBLOCK %s\n\n' % keyword)

        keyword = 'SPECIES_POT'
        fd.write('%%BLOCK %s\n' % keyword)
        for sp in sorted(self.pseudos):
            fd.write('    %s "%s"\n' % (sp[0], sp[1]))
        fd.write('%%ENDBLOCK %s\n\n' % keyword)

        keyword = 'SPECIES_ATOMIC_SET'
        fd.write('%%BLOCK %s\n' % keyword)
        for sp in sorted(self.solvers):
            fd.write('    %s "%s"\n' % (sp[0], sp[1]))
        fd.write('%%ENDBLOCK %s\n\n' % keyword)

        if ((self.parameters['ngwf_radius_cond'] > 0) or
            len(self.solvers_cond)==len(self.species)):
            keyword = 'SPECIES_ATOMIC_SET_COND'
            fd.write('%%BLOCK %s\n' % keyword)
            for sp in sorted(self.solvers_cond):
                fd.write('    %s "%s"\n' % (sp[0], sp[1]))
            fd.write('%%ENDBLOCK %s\n\n' % keyword)

        if self.core_wfs:
            keyword = 'SPECIES_CORE_WF'
            fd.write('%%BLOCK %s\n' % keyword)
            for sp in sorted(self.core_wfs):
                fd.write('    %s "%s"\n' % (sp[0], sp[1]))
            fd.write('%%ENDBLOCK %s\n\n' % keyword)

        if 'bsunfld_calculate' in self.parameters:
            if 'species_bsunfld_groups' not in self.parameters:
                self.parameters['species_bsunfld_groups'] = self.atoms.get_chemical_symbols()

        # Loop over parameters entries in alphabetal order, outputting
        # them as keywords or blocks as appropriate
        for p, param in sorted(parameters.items()):
            if param is not None and \
                    p.lower() not in self._dummy_parameters:
                if p.lower() in self._block_parameters:
                    keyword = p.upper()
                    fd.write('\n%%BLOCK %s\n' % keyword)
                    if p.lower() in self._path_parameters:
                        self.write_kpt_path(fd, param)
                    elif p.lower() in self._group_parameters:
                        self.write_groups(fd, param)
                    else:
                        fd.write('%s\n' % str(param))
                    fd.write('%%ENDBLOCK %s\n\n' % keyword)
                else:
                    fd.write('%s : %s\n' % (p, param))
            if p.upper() == 'XC':
                # Onetep calls XC something else...
                fd.write('xc_functional : %s\n' % param)
        fd.close()

    def write_kpt_path(self, fd, path):
        """Writes a k-point path to a ONETEP input file"""
        for kpt in array(path):
            fd.write('    %8.6f %8.6f %8.6f\n' % (kpt[0], kpt[1], kpt[2]))

    def write_groups(self, fd, groups):
        """Writes multiple groups of atom labels to a ONETEP input file"""
        for grp in groups:
            fd.write(" ".join(map(str, grp)))
            fd.write('\n')

    def __repr__(self):
        """Returns generic, fast to capture representation of
        ONETEP settings along with atoms object.
        """
        expr = ''
        expr += '-----------------Atoms--------------------\n'
        if self.atoms is not None:
            expr += str('%20s\n' % self.atoms)
        else:
            expr += 'None\n'

        expr += '\n-----------------Species---------------------\n'
        expr += str(self.species)
        expr += '\n-----------------Pseudos---------------------\n'
        expr += str(self.pseudos)
        expr += '\n-----------------Options------------\n'
        for key in self.parameters:
            expr += '%20s : %s\n' % (key, self.parameters[key])

        return expr

    def set_label(self, label):
        """The label is part of each seed, which in turn is a prefix
        in each ONETEP related file.
        """
        self.label = label
        self.prefix = label
