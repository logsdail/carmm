
# lammps.py (2011/03/29)
# An ASE calculator for the LAMMPS classical MD code available from
#       http://lammps.sandia.gov/
# The environment variable ASE_LAMMPSRUN_COMMAND must be defined to point to the
# LAMMPS binary.
#
# Copyright (C) 2009 - 2011 Joerg Meyer, joerg.meyer@ch.tum.de
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
# USA or see <http://www.gnu.org/licenses/>.


import os
import shutil
import shlex
from subprocess import Popen, PIPE
from threading import Thread
from re import compile as re_compile, IGNORECASE
from tempfile import mkdtemp, NamedTemporaryFile, mktemp as uns_mktemp
import inspect
import warnings
import numpy as np

from ase import Atoms
from ase.parallel import paropen
from ase.calculators.calculator import Calculator
from ase.calculators.calculator import all_changes
from ase.utils import basestring as asestring
from ase.data import chemical_symbols
from ase.data import atomic_masses
from ase.io.lammpsdata import write_lammps_data
from ase.io.lammpsrun import read_lammps_dump
from ase.calculators.lammps import Prism
from ase.calculators.lammps import write_lammps_in
from ase.calculators.lammps import CALCULATION_END_MARK
from ase.calculators.lammps import convert

__all__ = ["LAMMPS"]


class LAMMPS(Calculator):
    """The LAMMPS calculators object

    files: list
        List of files typically containing relevant potentials for the calculation
    parameters: dict
        Dictionary of settings to be passed into the input file for calculation.
    specorder: list
        Within LAAMPS, atoms are identified by an integer value starting from 1.
        This variable allows the user to define the order of the indices assigned to the
        atoms in the calculation, with the default if not given being alphabetical
    keep_tmp_files: bool
        Retain any temporary files created. Mostly useful for debugging.
    tmp_dir: str
        path/dirname (default None -> create automatically).
        Explicitly control where the calculator object should create
        its files. Using this option implies 'keep_tmp_files'
    no_data_file: bool
        Controls whether an explicit data file will be used for feeding
        atom coordinates into lammps. Enable it to lessen the pressure on
        the (tmp) file system. THIS OPTION MIGHT BE UNRELIABLE FOR CERTAIN
        CORNER CASES (however, if it fails, you will notice...).
    keep_alive: bool
        When using LAMMPS as a spawned subprocess, keep the subprocess
        alive (but idling when unused) along with the calculator object.
    always_triclinic: bool
        Force use of a triclinic cell in LAMMPS, even if the cell is
        a perfect parallelepiped.

        **Example**

Provided that the respective potential file is in the working directory, one
can simply run (note that LAMMPS needs to be compiled to work with EAM
potentials)

::

    from ase import Atom, Atoms
    from ase.build import bulk
    from ase.calculators.lammpsrun import LAMMPS

    parameters = {'pair_style': 'eam/alloy',
                  'pair_coeff': ['* * NiAlH_jea.eam.alloy H Ni']}

    files = ['NiAlH_jea.eam.alloy']

    Ni = bulk('Ni', cubic=True)
    H = Atom('H', position=Ni.cell.diagonal()/2)
    NiH = Ni + H

    lammps = LAMMPS(parameters=parameters, files=files)

    NiH.set_calculator(lammps)
    print("Energy ", NiH.get_potential_energy())

(Remember you also need to set the environment variable ``$ASE_LAMMPSRUN_COMMAND``)

    """

    name = "lammpsrun"
    implemented_properties = ["energy", "forces", "stress", "energies"]

    # parameters to choose options in LAMMPSRUN
    ase_parameters = dict(
        specorder=None,
        always_triclinic=False,
        keep_alive=True,
        keep_tmp_files=False,
        no_data_file=False,
        tmp_dir=None,
        files=[],  # usually contains potential parameters
        verbose=False,
        write_velocities=False,
        binary_dump=True,  # bool - use binary dump files (full
                           # precision but long long ids are casted to
                           # double)
        lammps_options="-echo log -screen none -log /dev/stdout",
        trajectory_out=None,  # file object, if is not None the
                              # trajectory will be saved in it
    )

    # parameters forwarded to LAMMPS
    lammps_parameters = dict(
        boundary=None,  # bounadry conditions styles
        units="metal",  # str - Which units used; some potentials
                        # require certain units
        atom_style="atomic",
        special_bonds=None,
        # potential informations
        pair_style="lj/cut 2.5",
        pair_coeff=["* * 1 1"],
        masses=None,
        pair_modify=None,
        # variables controlling the output
        thermo_args=[
            "step", "temp", "press", "cpu", "pxx", "pyy", "pzz",
            "pxy", "pxz", "pyz", "ke", "pe", "etotal", "vol", "lx",
            "ly", "lz", "atoms", ],
        dump_properties=["id", "type", "x", "y", "z", "vx", "vy",
                         "vz", "fx", "fy", "fz", ],
        dump_period=1,  # period of system snapshot saving (in MD steps)
    )

    default_parameters = dict(ase_parameters, **lammps_parameters)

    # legacy parameter persist, when the 'parameters' dictinary is manipulated
    # from the outside.  All others are rested to the default value
    legacy_parameters = [
        "specorder",
        "dump_period",
        "always_triclinic",
        "keep_alive",
        "keep_tmp_files",
        "tmp_dir",
        "parameters",
        "no_data_file",
        "files",
        "write_velocities",
        "trajectory_out",
    ]

    legacy_parameters_map = {"_custom_thermo_args": "thermo_args"}

    legacy_warn_string = "You are using an "
    legacy_warn_string += "old syntax to set '{}'.\n"
    legacy_warn_string += "Please use {}.set().".format(name.upper())

    def __init__(self, label="lammps", **kwargs):
        # "Parameters" used to be the dictionary with all parameters forwarded
        # to lammps.  This clashes with the implementation in Calculator to
        # reload an old one. Trying to catch both cases to not break old
        # scripts.
        if "parameters" in kwargs:
            old_parameters = kwargs["parameters"]
            if isinstance(old_parameters, dict):
                warnings.warn(self.legacy_warn_string.format("parameters"))
                del kwargs["parameters"]
        else:
            old_parameters = None

        Calculator.__init__(self, label=label, **kwargs)

        if old_parameters and isinstance(old_parameters, dict):
            self.set(**old_parameters)

        self.prism = None
        self.calls = 0
        self.forces = None
        # thermo_content contains data "written by" thermo_style.
        # It is a list of dictionaries, each dict (one for each line
        # printed by thermo_style) contains a mapping between each
        # custom_thermo_args-argument and the corresponding
        # value as printed by lammps. thermo_content will be
        # re-populated by the read_log method.
        self.thermo_content = []

        if self.parameters.tmp_dir is not None:
            # If tmp_dir is pointing somewhere, don't remove stuff!
            self.parameters.keep_tmp_files = True
        self._lmp_handle = None  # To handle the lmp process

        if self.parameters.tmp_dir is None:
            self.parameters.tmp_dir = mkdtemp(prefix="LAMMPS-")
        else:
            self.parameters.tmp_dir = os.path.realpath(self.parameters.tmp_dir)
            if not os.path.isdir(self.parameters.tmp_dir):
                os.mkdir(self.parameters.tmp_dir, 0o755)

        for f in self.parameters.files:
            shutil.copy(
                f, os.path.join(self.parameters.tmp_dir, os.path.basename(f))
            )

    def get_lammps_command(self):
        cmd = self.parameters.get('command')
        if cmd is None:
            envvar = 'ASE_{}_COMMAND'.format(self.name.upper())
            cmd = os.environ.get(envvar)

        if cmd is None:
            cmd = 'lammps'

        opts = self.parameters.get('lammps_options')

        if opts is not None:
            cmd = '{} {}'.format(cmd, opts)

        return cmd

    def __setattr__(self, key, value):
        """Catch attribute sets to emulate legacy behavior.

        Old LAMMPSRUN allows to just override the parameters
        dictionary. "Modern" ase calculators can assume that default
        parameters are always set, overrides of the
        'parameters'-dictionary have to be caught and the default
        parameters need to be added first.  A check refuses to set
        calculator attributes if they are unknown and set outside the
        '__init__' functions.
        """
        # !TODO: remove and break somebody's code (e.g. the test examples)
        if (
                key == "parameters"
                and value is not None
                and self.parameters is not None
        ):
            temp_dict = self.get_default_parameters()
            if self.parameters:
                for l_key in self.legacy_parameters:
                    try:
                        temp_dict[l_key] = self.parameters[l_key]
                    except KeyError:
                        pass
            temp_dict.update(value)
            value = temp_dict
        if key in self.legacy_parameters and key != "parameters":
            warnings.warn(self.legacy_warn_string.format(key))
            self.set(**{key: value})
        elif key in self.legacy_parameters_map:
            warnings.warn(
                self.legacy_warn_string.format(
                    "{} for {}".format(self.legacy_parameters_map[key], key)
                )
            )
            self.set(**{self.legacy_parameters_map[key]: value})
        # Catch setting none-default attributes
        # one test was assigning an useless Attribute, but it still worked
        # because the assigned object was before manipulation already handed
        # over to the calculator (10/2018)
        elif hasattr(self, key) or inspect.stack()[1][3] == "__init__":
            Calculator.__setattr__(self, key, value)
        else:
            raise AttributeError("Setting unknown Attribute '{}'".format(key))

    def __getattr__(self, key):
        """Corresponding getattribute-function to emulate legacy behavior.
        """
        if key in self.legacy_parameters and key != "parameters":
            return self.parameters[key]
        if key in self.legacy_parameters_map:
            return self.parameters[self.legacy_parameters_map[key]]
        return object.__getattribute__(self, key)

    def clean(self, force=False):

        self._lmp_end()

        if not self.parameters.keep_tmp_files or force:
            shutil.rmtree(self.parameters.tmp_dir)

    def check_state(self, atoms, tol=1.0e-10):
        # Transforming the unit cell to conform to LAMMPS' convention for
        # orientation (c.f. https://lammps.sandia.gov/doc/Howto_triclinic.html)
        # results in some precision loss, so we use bit larger tolerance than
        # machine precision here.  Note that there can also be precision loss
        # related to how many significant digits are specified for things in
        # the LAMMPS input file.
        return Calculator.check_state(self, atoms, tol)

    def calculate(self, atoms=None, properties=None, system_changes=None):
        if properties is None:
            properties = self.implemented_properties
        if system_changes is None:
            system_changes = all_changes
        Calculator.calculate(self, atoms, properties, system_changes)
        self.run()

    def _lmp_alive(self):
        # Return True if this calculator is currently handling a running
        # lammps process
        return self._lmp_handle and not isinstance(
            self._lmp_handle.poll(), int
        )

    def _lmp_end(self):
        # Close lammps input and wait for lammps to end. Return process
        # return value
        if self._lmp_alive():
            self._lmp_handle.stdin.close()
            # !TODO: handle lammps error codes
            # return self._lmp_handle.wait()

    def set_missing_parameters(self):
        """Verify that all necessary variables are set.
        """
        symbols = self.atoms.get_chemical_symbols()
        # If unspecified default to atom types in alphabetic order
        if not self.parameters.specorder:
            self.parameters.specorder = sorted(set(symbols))

        # !TODO: handle cases were setting masses actual lead to errors
        if not self.parameters.masses:
            self.parameters.masses = []
            for type_id, specie in enumerate(self.parameters.specorder):
                mass = atomic_masses[chemical_symbols.index(specie)]
                self.parameters.masses += [
                    "{0:d} {1:f}".format(type_id + 1, mass)
                ]

        # set boundary condtions
        if not self.parameters.boundary:
            b_str = " ".join(["fp"[int(x)] for x in self.atoms.get_pbc()])
            self.parameters.boundary = b_str

    def run(self, set_atoms=False):
        # !TODO: split this function
        """Method which explicitly runs LAMMPS."""
        pbc = self.atoms.get_pbc()
        if all(pbc):
            cell = self.atoms.get_cell()
        elif not any(pbc):
            # large enough cell for non-periodic calculation -
            # LAMMPS shrink-wraps automatically via input command
            #       "periodic s s s"
            # below
            cell = 2 * np.max(np.abs(self.atoms.get_positions())) * np.eye(3)
        else:
            warnings.warn(
                "semi-periodic ASE cell detected - translation "
                + "to proper LAMMPS input cell might fail"
            )
            cell = self.atoms.get_cell()
        self.prism = Prism(cell)

        self.set_missing_parameters()
        self.calls += 1

        # change into subdirectory for LAMMPS calculations
        cwd = os.getcwd()
        os.chdir(self.parameters.tmp_dir)

        # setup file names for LAMMPS calculation
        label = "{0}{1:>06}".format(self.label, self.calls)
        lammps_in = uns_mktemp(
            prefix="in_" + label, dir=self.parameters.tmp_dir
        )
        lammps_log = uns_mktemp(
            prefix="log_" + label, dir=self.parameters.tmp_dir
        )
        lammps_trj_fd = NamedTemporaryFile(
            prefix="trj_" + label,
            suffix=(".bin" if self.parameters.binary_dump else ""),
            dir=self.parameters.tmp_dir,
            delete=(not self.parameters.keep_tmp_files),
        )
        lammps_trj = lammps_trj_fd.name
        if self.parameters.no_data_file:
            lammps_data = None
        else:
            lammps_data_fd = NamedTemporaryFile(
                prefix="data_" + label,
                dir=self.parameters.tmp_dir,
                delete=(not self.parameters.keep_tmp_files),
                mode='w',
                encoding='ascii'
            )
            write_lammps_data(
                lammps_data_fd,
                self.atoms,
                specorder=self.parameters.specorder,
                force_skew=self.parameters.always_triclinic,
                velocities=self.parameters.write_velocities,
                prismobj=self.prism,
                units=self.parameters.units,
                atom_style=self.parameters.atom_style
            )
            lammps_data = lammps_data_fd.name
            lammps_data_fd.flush()

        # see to it that LAMMPS is started
        if not self._lmp_alive():
            command = self.get_lammps_command()
            # Attempt to (re)start lammps
            self._lmp_handle = Popen(
                shlex.split(command, posix=(os.name == "posix")),
                stdin=PIPE,
                stdout=PIPE,
            )
        lmp_handle = self._lmp_handle

        # Create thread reading lammps stdout (for reference, if requested,
        # also create lammps_log, although it is never used)
        if self.parameters.keep_tmp_files:
            lammps_log_fd = open(lammps_log, "wb")
            fd = SpecialTee(lmp_handle.stdout, lammps_log_fd)
        else:
            fd = lmp_handle.stdout
        thr_read_log = Thread(target=self.read_lammps_log, args=(fd,))
        thr_read_log.start()

        # write LAMMPS input (for reference, also create the file lammps_in,
        # although it is never used)
        if self.parameters.keep_tmp_files:
            lammps_in_fd = open(lammps_in, "wb")
            fd = SpecialTee(lmp_handle.stdin, lammps_in_fd)
        else:
            fd = lmp_handle.stdin
        write_lammps_in(
            lammps_in=fd,
            parameters=self.parameters,
            atoms=self.atoms,
            prismobj=self.prism,
            lammps_trj=lammps_trj,
            lammps_data=lammps_data,
        )

        if self.parameters.keep_tmp_files:
            lammps_in_fd.close()

        # Wait for log output to be read (i.e., for LAMMPS to finish)
        # and close the log file if there is one
        thr_read_log.join()
        if self.parameters.keep_tmp_files:
            lammps_log_fd.close()

        if not self.parameters.keep_alive:
            self._lmp_end()

        exitcode = lmp_handle.poll()
        if exitcode and exitcode != 0:
            cwd = os.getcwd()
            raise RuntimeError(
                "LAMMPS exited in {} with exit code: {}."
                "".format(cwd, exitcode)
            )

        # A few sanity checks
        if len(self.thermo_content) == 0:
            raise RuntimeError("Failed to retrieve any thermo_style-output")
        if int(self.thermo_content[-1]["atoms"]) != len(self.atoms):
            # This obviously shouldn't happen, but if prism.fold_...() fails,
            # it could
            raise RuntimeError("Atoms have gone missing")

        trj_atoms = read_lammps_dump(
            infileobj=lammps_trj,
            order=False,
            index=-1,
            prismobj=self.prism,
            specorder=self.parameters.specorder,
        )

        if set_atoms:
            self.atoms = trj_atoms.copy()

        self.forces = trj_atoms.get_forces()
        # !TODO: trj_atoms is only the last snapshot of the system; Is it
        #        desireable to save also the inbetween steps?
        if self.parameters.trajectory_out is not None:
            # !TODO: is it advisable to create here temporary atoms-objects
            self.trajectory_out.write(trj_atoms)

        tc = self.thermo_content[-1]
        self.results["energy"] = convert(
            tc["pe"], "energy", self.parameters["units"], "ASE"
        )
        self.results["free_energy"] = self.results["energy"]
        self.results["forces"] = self.forces.copy()
        stress = np.array(
            [-tc[i] for i in ("pxx", "pyy", "pzz", "pyz", "pxz", "pxy")]
        )

        # We need to apply the Lammps rotation stuff to the stress:
        xx, yy, zz, yz, xz, xy = stress
        stress_tensor = np.array([[xx, xy, xz],
                                  [xy, yy, yz],
                                  [xz, yz, zz]])
        R = self.prism.rot_mat
        stress_atoms = np.dot(R, stress_tensor)
        stress_atoms = np.dot(stress_atoms, R.T)
        stress_atoms = stress_atoms[[0, 1, 2, 1, 0, 0],
                                    [0, 1, 2, 2, 2, 1]]
        stress = stress_atoms

        self.results["stress"] = convert(
            stress, "pressure", self.parameters["units"], "ASE"
        )

        lammps_trj_fd.close()
        if not self.parameters.no_data_file:
            lammps_data_fd.close()

        os.chdir(cwd)

    def read_lammps_log(self, lammps_log=None):
        # !TODO: somehow communicate 'thermo_content' explicitly
        """Method which reads a LAMMPS output log file."""

        if lammps_log is None:
            lammps_log = self.label + ".log"

        if isinstance(lammps_log, asestring):
            fileobj = paropen(lammps_log, "wb")
            close_log_file = True
        else:
            # Expect lammps_in to be a file-like object
            fileobj = lammps_log
            close_log_file = False

        # read_log depends on that the first (three) thermo_style custom args
        # can be capitilized and matched against the log output. I.e.
        # don't use e.g. 'ke' or 'cpu' which are labeled KinEng and CPU.
        _custom_thermo_mark = " ".join(
            [x.capitalize() for x in self.parameters.thermo_args[0:3]]
        )

        # !TODO: regex-magic necessary?
        # Match something which can be converted to a float
        f_re = r"([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?|nan|inf))"
        n_args = len(self.parameters["thermo_args"])
        # Create a re matching exactly N white space separated floatish things
        _custom_thermo_re = re_compile(
            r"^\s*" + r"\s+".join([f_re] * n_args) + r"\s*$", flags=IGNORECASE
        )

        thermo_content = []
        line = fileobj.readline().decode("utf-8")
        while line and line.strip() != CALCULATION_END_MARK:
            # check error
            if 'ERROR:' in line:
                if close_log_file:
                    fileobj.close()
                raise RuntimeError('LAMMPS exits with error message: {}'.format(line))

            # get thermo output
            if line.startswith(_custom_thermo_mark):
                bool_match = True
                while bool_match:
                    line = fileobj.readline().decode("utf-8")
                    bool_match = _custom_thermo_re.match(line)
                    if bool_match:
                        # create a dictionary between each of the
                        # thermo_style args and it's corresponding value
                        thermo_content.append(
                            dict(
                                zip(
                                    self.parameters.thermo_args,
                                    map(float, bool_match.groups()),
                                )
                            )
                        )
            else:
                line = fileobj.readline().decode("utf-8")

        if close_log_file:
            fileobj.close()

        self.thermo_content = thermo_content


class SpecialTee:
    """A special purpose, with limited applicability, tee-like thing.

    A subset of stuff read from, or written to, orig_fd,
    is also written to out_fd.
    It is used by the lammps calculator for creating file-logs of stuff
    read from, or written to, stdin and stdout, respectively.
    """

    def __init__(self, orig_fd, out_fd):
        self._orig_fd = orig_fd
        self._out_fd = out_fd
        self.name = orig_fd.name

    def write(self, data):
        self._orig_fd.write(data)
        self._out_fd.write(data)
        self.flush()

    def read(self, *args, **kwargs):
        data = self._orig_fd.read(*args, **kwargs)
        self._out_fd.write(data)
        return data

    def readline(self, *args, **kwargs):
        data = self._orig_fd.readline(*args, **kwargs)
        self._out_fd.write(data)
        return data

    def readlines(self, *args, **kwargs):
        data = self._orig_fd.readlines(*args, **kwargs)
        self._out_fd.write("".join(data))
        return data

    def flush(self):
        self._orig_fd.flush()
        self._out_fd.flush()


if __name__ == "__main__":
    pair_style = "eam"
    Pd_eam_file = "Pd_u3.eam"
    pair_coeff = ["* * " + Pd_eam_file]
    parameters = {"pair_style": pair_style, "pair_coeff": pair_coeff}
    my_files = [Pd_eam_file]
    calc = LAMMPS(parameters=parameters, files=my_files)
    a0 = 3.93
    b0 = a0 / 2.0

    bulk = Atoms(
        ["Pd"] * 4,
        positions=[(0, 0, 0), (b0, b0, 0), (b0, 0, b0), (0, b0, b0)],
        cell=[a0] * 3,
        pbc=True,
    )
    # test get_forces
    print("forces for a = {0}".format(a0))
    print(calc.get_forces(bulk))
    # single points for various lattice constants
    bulk.set_calculator(calc)
    for i in range(-5, 5, 1):
        a = a0 * (1 + i / 100.0)
        bulk.set_cell([a] * 3)
        print("a : {0} , total energy : {1}".format(
            a, bulk.get_potential_energy()))

    calc.clean()
