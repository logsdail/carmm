# run using python3 <filename>.py
from ase.build import fcc111, fcc100, add_adsorbate
from ase.visualize import view
from ase.io import write
from math import sqrt

# Edit these
atomic_species='Au'
lattice_parameter=2.939
unit_cell_depth=4
unit_cell_width=4
slab_depth=4
vacuum_region_size=10.0

# Create surface
slab = fcc111(atomic_species, a=lattice_parameter*sqrt(2), size=(unit_cell_width,unit_cell_depth,slab_depth))
#slab = fcc100(atomic_species, a=lattice_parameter*sqrt(2), size=(unit_cell_width,unit_cell_depth,slab_depth))

# Enable to add H at an ontop position
#add_adsorbate(slab, 'H', 1.5, 'ontop')

# Add vacuum
slab.center(vacuum=vacuum_region_size, axis=2)

# Allows you to visualise the slab
view(slab)
# Writes surface to file
#write('slab_geometry.in',slab,format='aims')
