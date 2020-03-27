def create_facets(bulk_model, layers=2, facets=[(1,1,1)], supercell =(1,1,1), vacuum=10):
	from ase.spacegroup import crystal
	from ase.io import write, read
	from ase.build import surface
	#input bulk geometry
	if isinstance(bulk_model, str) is True:
		bulk_model = read(bulk_model)
	#loop for all surfaces
	for face in facets:
		slab = surface(bulk_model,(face),layers)		#cut the slab in face, last number is the number of layers
		slab.repeat(supercell)					#repeat the slab in x,y,z direction and create a supercell
		slab.center(vacuum=vacuum, axis=2)			#introduce vaccum around the supercell in y direction

		#write into file, name based on facet type
		write(str(face[0])+str(face[1])+str(face[2])+'.in', slab)
