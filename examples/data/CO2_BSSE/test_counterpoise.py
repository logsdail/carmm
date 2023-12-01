from ase.io import read
from carmm.analyse.counterpoise_onepot import counterpoise_calc
from carmm.run.aims_calculator import get_aims_calculator
from carmm.run.aims_path import set_aims_command
import subprocess

CO = read('C_monoxide_pbe.traj')
set_aims_command(hpc='hawk', basis_set='light', defaults=2020)
toy_calc = get_aims_calculator(dimensions=0, xc='pbe')
toy_calc.set(xc='pbe', spin='collinear', default_initial_moment=0.5, relativistic='atomic_zora scalar')
toy_calc.set(species_dir='.')
cp_index = counterpoise_calc(CO, a_id=[1], b_id=[0], symbol_not_index=False, fhi_calc=toy_calc, dry_run=True)
cp_symbol = counterpoise_calc(CO, a_id=['C'], b_id=['O'], symbol_not_index=True, fhi_calc=toy_calc, dry_run=True)
assert cp_index == -8.707358006176946e-05
assert cp_symbol == -8.707358006176946e-05

f = open("geometry.in", 'r')
lines = f.readlines()
assert lines[6] == "empty -0.0000000000000000 0.0000000000000000 -0.6536947973321450 C\n"
subprocess.check_call(['rm', 'control.in'])
subprocess.check_call(['rm', 'geometry.in'])
subprocess.check_call(['rm', 'parameters.ase'])

