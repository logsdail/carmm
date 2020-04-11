from pathlib import Path
import numpy as np
from ase.build import bulk
from ase.io.aims import read_aims as read

parent = Path(__file__).parent
format = "aims"

atoms = bulk("Si")
atoms.positions[0, 0] -= 0.01

file = "geometry.in"


# check cartesian
def test_cartesian(atoms=atoms):
    """write cartesian coords and check if structure was preserved"""
    atoms.write(file, format=format)

    new_atoms = read((file))

    assert np.allclose(atoms.positions, new_atoms.positions)

    Path(file).unlink()


# check scaled
def test_scaled(atoms=atoms):
    """write fractional coords and check if structure was preserved"""
    atoms.write(file, format=format, scaled=True, wrap=False)

    new_atoms = read(file)

    assert np.allclose(atoms.positions, new_atoms.positions), (
        atoms.positions,
        new_atoms.positions,
    )

    Path(file).unlink()


# this should fail
def test_scaled_wrapped(atoms=atoms):
    """write fractional coords and check if structure was preserved"""
    atoms.write(file, format=format, scaled=True, wrap=True)

    new_atoms = read(file)

    try:
        assert np.allclose(atoms.positions, new_atoms.positions), (
            atoms.positions,
            new_atoms.positions,
        )
    except AssertionError:
        atoms.wrap()
        assert np.allclose(atoms.positions, new_atoms.positions), (
            atoms.positions,
            new_atoms.positions,
        )

    Path(file).unlink()


if __name__ == "__main__":
    test_cartesian()
    test_scaled()
    test_scaled_wrapped()
