#TODO: How about we change facets to "highest_index_facet", and it just generates all possible facets up to this number?
#      I'm also going to change this so it doesn't write by default, as I'd prefer it to return an object
def generate(bulk_model, layers=2, facets=[(1,1,1)], supercell =(1,1,1), vacuum=10, save=False):
	'''
	Function to create variety of different facets

	Parameters:

	bulk_model: Atoms object or string
		*Description*
	layers: int
		*Description*
	facets: list of array of ints
		*Description - does this need to be a list?*
	supercell: array of ints
		*Description - can this just be an int?*
	vacuum: float
		*Description*

	'''

	from ase.build import surface

	#input bulk geometry
	#if string, check if atomic label, in which case we don't read from file
	if isinstance(bulk_model, str) and len(bulk_model) > 2:
		from ase.io import read
		bulk_model = read(bulk_model)

	slabs = []
	#loop for all surfaces
	for face in facets:
		# cut the slab in face, last number is the number of layers
		slab = surface(bulk_model,(face),layers)
		# repeat the slab in x,y,z direction and create a supercell
		slab.repeat(supercell)
		# introduce vaccum around the supercell in y direction - do you mean z?
		slab.center(vacuum=vacuum, axis=2)
		# Add model to the list of slabs
		slabs.append(slab)

	# Now writes all files at thend
	if save:
		_save(facets, slabs)

	return facets, slabs

def _save(facets, models):
	'''
	Method to write all structures to file; hidden from user

	Parameters:

	facets: list of array of ints
		*Description - this now needs to be a list*
	slabs: list of atom objects
		*Description*

	'''

	from ase.io import write

	for i in range(len(facets)):
		face = facets[i]
		write(str(face[0])+str(face[1])+str(face[2])+'.in', models[i])
