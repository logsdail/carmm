def create_facets(bulk_model, layers=2, facets=[(1,1,1)], supercell =(1,1,1), vacuum=10):
	'''
	Function to create variety of different facets

	Parameters:

	bulk_model: Atoms object or string
		*Description*
	layers: int
		*Description*
	facets: list of array of ints
		*Description - does this need to be a list?
	supercell: array of ints
		*Description - can this just be an int?*
	vacuum: float
		*Description

	'''

	from ase.io import write
	from ase.build import surface

	#input bulk geometry
	if isinstance(bulk_model, str):
		from ase.io import read
		bulk_model = read(bulk_model)

	#loop for all surfaces
	for face in facets:
		# cut the slab in face, last number is the number of layers
		slab = surface(bulk_model,(face),layers)
		# repeat the slab in x,y,z direction and create a supercell
		slab.repeat(supercell)
		# introduce vaccum around the supercell in y direction - do you mean z?
		slab.center(vacuum=vacuum, axis=2)

		#TODO: Does this need to include the write function? Could it just return the results to be saved elsewhere?
		#      How about just returning the list of facets at the end, and leaving the user to use the write function?
		#write into file, name based on facet type
		write(str(face[0])+str(face[1])+str(face[2])+'.in', slab)
