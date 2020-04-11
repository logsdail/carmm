"""Structure optimization. """

import sys
import pickle
import time
from math import sqrt
from os.path import isfile

from ase.calculators.calculator import PropertyNotImplementedError
from ase.parallel import world, barrier
from ase.io.trajectory import Trajectory
import collections


class Dynamics:
    """Base-class for all MD and structure optimization classes."""

    def __init__(
        self, atoms, logfile, trajectory, append_trajectory=False, master=None
    ):
        """Dynamics object.

        Parameters:

        atoms: Atoms object
            The Atoms object to operate on.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

        trajectory: Trajectory object or str
            Attach trajectory object.  If *trajectory* is a string a
            Trajectory will be constructed.  Use *None* for no
            trajectory.

        append_trajectory: boolean
            Defaults to False, which causes the trajectory file to be
            overwriten each time the dynamics is restarted from scratch.
            If True, the new structures are appended to the trajectory
            file instead.

        master: boolean
            Defaults to None, which causes only rank 0 to save files.  If
            set to true,  this rank will save files.
        """

        self.atoms = atoms
        if master is None:
            master = world.rank == 0
        if not master:
            logfile = None
        elif isinstance(logfile, str):
            if logfile == "-":
                logfile = sys.stdout
            else:
                logfile = open(logfile, "a")
        self.logfile = logfile

        self.observers = []
        self.nsteps = 0
        # maximum number of steps placeholder with maxint
        self.max_steps = 100000000

        if trajectory is not None:
            if isinstance(trajectory, str):
                mode = "a" if append_trajectory else "w"
                trajectory = Trajectory(
                    trajectory, mode=mode, atoms=atoms, master=master
                )
            self.attach(trajectory)

    def get_number_of_steps(self):
        return self.nsteps

    def insert_observer(
        self, function, position=0, interval=1, *args, **kwargs
    ):
        """Insert an observer."""
        if not isinstance(function, collections.Callable):
            function = function.write
        self.observers.insert(position, (function, interval, args, kwargs))

    def attach(self, function, interval=1, *args, **kwargs):
        """Attach callback function.

        If *interval > 0*, at every *interval* steps, call *function* with
        arguments *args* and keyword arguments *kwargs*.

        If *interval <= 0*, after step *interval*, call *function* with
        arguments *args* and keyword arguments *kwargs*.  This is
        currently zero indexed."""

        if hasattr(function, "set_description"):
            d = self.todict()
            d.update(interval=interval)
            function.set_description(d)
        if not hasattr(function, "__call__"):
            function = function.write
        self.observers.append((function, interval, args, kwargs))

    def call_observers(self):
        for function, interval, args, kwargs in self.observers:
            call = False
            # Call every interval iterations
            if interval > 0:
                if (self.nsteps % interval) == 0:
                    call = True
            # Call only on iteration interval
            elif interval <= 0:
                if self.nsteps == abs(interval):
                    call = True
            if call:
                function(*args, **kwargs)

    def irun(self):
        """Run dynamics algorithm as generator. This allows, e.g.,
        to easily run two optimizers or MD thermostats at the same time.

        Examples:
        >>> opt1 = BFGS(atoms)
        >>> opt2 = BFGS(StrainFilter(atoms)).irun()
        >>> for _ in opt2:
        >>>     opt1.run()
        """

        # compute inital structure and log the first step
        self.atoms.get_forces()

        # yield the first time to inspect before logging
        yield False

        if self.nsteps == 0:
            self.log()
            self.call_observers()

        # run the algorithm until converged or max_steps reached
        while not self.converged() and self.nsteps < self.max_steps:

            # compute the next step
            self.step()
            self.nsteps += 1

            # let the user inspect the step and change things before logging
            # and predicting the next step
            yield False

            # log the step
            self.log()
            self.call_observers()

        # finally check if algorithm was converged
        yield self.converged()

    def run(self):
        """Run dynamics algorithm.

        This method will return when the forces on all individual
        atoms are less than *fmax* or when the number of steps exceeds
        *steps*."""

        for converged in Dynamics.irun(self):
            pass
        return converged

    def converged(self, *args):
        """" a dummy function as placeholder for a real criterion, e.g. in
        Optimizer """
        return False

    def log(self, *args):
        """ a dummy function as placeholder for a real logger, e.g. in
        Optimizer """
        return True

    def step(self):
        """this needs to be implemented by subclasses"""
        raise RuntimeError("step not implemented.")


class Optimizer(Dynamics):
    """Base-class for all structure optimization classes."""

    def __init__(
        self,
        atoms,
        restart,
        logfile,
        trajectory,
        master=None,
        append_trajectory=False,
        force_consistent=False,
    ):
        """Structure optimizer object.

        Parameters:

        atoms: Atoms object
            The Atoms object to relax.

        restart: str
            Filename for restart file.  Default value is *None*.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

        trajectory: Trajectory object or str
            Attach trajectory object.  If *trajectory* is a string a
            Trajectory will be constructed.  Use *None* for no
            trajectory.

        master: boolean
            Defaults to None, which causes only rank 0 to save files.  If
            set to true,  this rank will save files.

        append_trajectory: boolean
            Appended to the trajectory file instead of overwriting it.

        force_consistent: boolean or None
            Use force-consistent energy calls (as opposed to the energy
            extrapolated to 0 K).  If force_consistent=None, uses
            force-consistent energies if available in the calculator, but
            falls back to force_consistent=False if not.
        """
        Dynamics.__init__(
            self,
            atoms,
            logfile,
            trajectory,
            append_trajectory=append_trajectory,
            master=master,
        )

        self.force_consistent = force_consistent
        if self.force_consistent is None:
            self.set_force_consistent()

        self.restart = restart

        # initialize attribute
        self.fmax = None

        if restart is None or not isfile(restart):
            self.initialize()
        else:
            self.read()
            barrier()

    def todict(self):
        description = {
            "type": "optimization",
            "optimizer": self.__class__.__name__,
        }
        return description

    def initialize(self):
        pass

    def irun(self, fmax=0.05, steps=None):
        """ call Dynamics.irun and keep track of fmax"""
        self.fmax = fmax
        if steps:
            self.max_steps = steps
        return Dynamics.irun(self)

    def run(self, fmax=0.05, steps=None):
        """ call Dynamics.run and keep track of fmax"""
        self.fmax = fmax
        if steps:
            self.max_steps = steps
        return Dynamics.run(self)

    def converged(self, forces=None):
        """Did the optimization converge?"""
        if forces is None:
            forces = self.atoms.get_forces()
        if hasattr(self.atoms, "get_curvature"):
            return (forces ** 2).sum(
                axis=1
            ).max() < self.fmax ** 2 and self.atoms.get_curvature() < 0.0
        return (forces ** 2).sum(axis=1).max() < self.fmax ** 2

    def log(self, forces=None):
        if forces is None:
            forces = self.atoms.get_forces()
        fmax = sqrt((forces ** 2).sum(axis=1).max())
        e = self.atoms.get_potential_energy(
            force_consistent=self.force_consistent
        )
        T = time.localtime()
        if self.logfile is not None:
            name = self.__class__.__name__
            if self.nsteps == 0:
                args = (" " * len(name), "Step", "Time", "Energy", "fmax")
                msg = "%s  %4s %8s %15s %12s\n" % args
                self.logfile.write(msg)

                if self.force_consistent:
                    msg = "*Force-consistent energies used in optimization.\n"
                    self.logfile.write(msg)

            ast = {1: "*", 0: ""}[self.force_consistent]
            args = (name, self.nsteps, T[3], T[4], T[5], e, ast, fmax)
            msg = "%s:  %3d %02d:%02d:%02d %15.6f%1s %12.4f\n" % args
            self.logfile.write(msg)

            self.logfile.flush()

    def dump(self, data):
        if world.rank == 0 and self.restart is not None:
            pickle.dump(data, open(self.restart, "wb"), protocol=2)

    def load(self):
        return pickle.load(open(self.restart, "rb"))

    def set_force_consistent(self):
        """Automatically sets force_consistent to True if force_consistent
        energies are supported by calculator; else False."""
        try:
            self.atoms.get_potential_energy(force_consistent=True)
        except PropertyNotImplementedError:
            self.force_consistent = False
        else:
            self.force_consistent = True
