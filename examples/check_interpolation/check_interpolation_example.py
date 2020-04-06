from software.analyse.neb_tools.check_geometry import switch_indices as Switch
from software.analyse.neb_tools.check_geometry import check_interpolation
from ase.io import read
index1 = 46
index2 = 47

changed = Switch('final.traj', index1, index2)
#initial = read("initial.traj")
#final = read("final.traj)"

#check_interpolation('initial.traj','final.traj',10)
check_interpolation('initial.traj', changed ,10)

#from ase.visualize import view
