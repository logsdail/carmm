"""Check that QE calculation can run."""

from ase.build import bulk
from ase.calculators.espresso import Espresso

# Default pseudos can go in ~/espresso/pseudo
# Get these from SSSP http://materialscloud.org/sssp/
#PSEUDO = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}

# I am updating this so we use a file which is provided by the
# Debian/Ubuntu package for QE.  Then it works out of the
# box, which is good for testing.  --askhl
PSEUDO = {'Si': 'Si.rel-pbe-rrkj.UPF'}
# /usr/share/espresso/pseudo/Si.rel-pbe-rrkj.UPF

# Don't forget to
# export ASE_ESPRESSO_COMMAND="mpirun -n 4 $HOME/Compile/q-e/bin/pw.x -in PREFIX.pwi > PREFIX.pwo"
# export ESPRESSO_PSEUDO="/path/to/pseudos"

def main():
    silicon = bulk('Si')
    calc = Espresso(pseudopotentials=PSEUDO,
                    ecutwfc=50.0)
    silicon.set_calculator(calc)
    silicon.get_potential_energy()

    assert calc.get_fermi_level() is not None
    assert calc.get_ibz_k_points() is not None
    assert calc.get_eigenvalues(spin=0, kpt=0) is not None
    assert calc.get_number_of_spins() is not None
    assert calc.get_k_point_weights() is not None

main()
