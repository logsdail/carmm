from ase.spacegroup import crystal
from ase.io import write,read
from ase.build import surface,supercells,cut
from ase.build.supercells import make_supercell
from ase.visualize import view

#input bulk geometry
bulk=read('bulk.in')

#list of faces
facets=[(0,0,1),(0,1,1), (1,1,1), (1,1,0), (1,0,0), (1,0,1), (0,1,0)]

#number of layers in z direction
layers = 1

#loop for all surfaces
for face in facets:
	slab = surface(bulk,(face),layers)		#cut the slab in face, last number is the number of layers
	slab.repeat((1,1,1))					#repeat the slab in x,y,z direction and create a supercell
	slab.center(vacuum=10, axis=2)			#introduce vaccum around the supercell in y direction
	write('print.in', slab)					#output the file 



## last line needs fixing to generate new name for each loop run. ##


