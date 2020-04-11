"""Test to ensure that md logging interval is respected (issue #304)."""

from ase.optimize import FIRE, BFGS
from ase.data import s22
from ase.calculators.tip3p import TIP3P
from ase.constraints import FixBondLengths
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase.io import Trajectory
import ase.units as u
import numpy as np
import os

md_cls_and_kwargs = [
    (VelocityVerlet, {}),
    (Langevin, {"temperature": 300 * u.kB, "friction": 0.02}),
]

def make_dimer(constraint=True):
    """Prepare atoms object for testing"""
    dimer = s22.create_s22_system("Water_dimer")

    calc = TIP3P(rc=9.0)
    dimer.calc = calc
    if constraint:
        dimer.constraints = FixBondLengths(
            [(3 * i + j, 3 * i + (j + 1) % 3) for i in range(2) for j in [0, 1, 2]]
            )
    return dimer

def fmax(forces):
    return np.sqrt((forces ** 2).sum(axis=1).max())


def test_opt(cls, atoms, logfile="opt.log", trajectory="opt.traj"):
    """run optimization and verify that log and trajectory have same number of steps"""

    # clean files to make sure the correct ones are tested
    for f in (logfile, trajectory):
        if os.path.exists(f):
            os.unlink(f)

    print("Testing", str(cls))
    opt = cls(atoms, logfile=logfile, trajectory=trajectory)

    # Run optimizer two times
    opt.run(0.2)
    opt.run(0.1)

    # Test number of lines in log file matches number of frames in trajectory
    with open(logfile, 'rt') as lf:
        lines = [l for l in lf]
    loglines = len(lines)
    print("Number of lines in log file:", loglines)

    with Trajectory(trajectory) as traj:
        trajframes = len(traj)
    print("Number of frames in trajectory:", trajframes)

    # There is a header line in the logfile
    assert loglines == trajframes + 1


def test_md(
    cls, atoms, kwargs, logfile="md.log", timestep=1 * u.fs, trajectory="md.traj", loginterval=1
):
    """ run MD for 10 steps and verify that log and trajectory have same number of steps"""

     # clean files to make sure the correct ones are tested
    for f in (logfile, trajectory):
        if os.path.exists(f):
            os.unlink(f)
            
    assert not atoms.constraints

    print("Testing", str(cls))
    md = cls(atoms, logfile=logfile, timestep=timestep, trajectory=trajectory, loginterval=loginterval, **kwargs)

    # run md two times
    md.run(steps=5)
    md.run(steps=5)

    # Test number of lines in log file matches number of frames in trajectory
    with open(logfile, 'rt') as lf:
        lines = [l for l in lf]
    loglines = len(lines)
    print("Number of lines in log file:", loglines)

    with Trajectory(trajectory) as traj:
        trajframes = len(traj)
    print("Number of frames in trajectory:", trajframes)

    # There is a header line in the logfile
    assert loglines == trajframes + 1

# test optimizer
for cls in (FIRE, BFGS):
    dimer = make_dimer()
    test_opt(cls, dimer)

# test md
for cls, kwargs in md_cls_and_kwargs:
    dimer = make_dimer(constraint=False)
    test_md(cls, dimer, kwargs=kwargs)

for cls, kwargs in md_cls_and_kwargs:
    dimer = make_dimer(constraint=False)
    test_md(cls, dimer, loginterval=2, kwargs=kwargs)
