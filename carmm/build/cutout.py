def cutout_sphere(atoms, centre, distance_cutoff=5.0):
    '''
    Returns a spherical cutout of a structure
    TODO: This doesn't work with periodic boundary conditions.
    It'd be nice to get that working and update the regression to test this.

    Parameters:

    atoms: Atoms object
        Input structure to cutout from
    centre: Integer
        Index of central atom in cutout
    distance_cutoff: Float
        Distance outside of which atoms are removed
    '''

    import numpy as np

    # Prevents unexpected editing of parent object in place
    # Now ensures returned object is different to incoming atoms
    atoms = atoms.copy()
    atoms_to_delete = []
    for i in range(len(atoms)):

        # get distances between atom of interest and others - then removes
        # all atoms beyond a certain radius

        #################### Edit atom tag ###############################
        distance_ab = np.linalg.norm((atoms.positions[centre] - atoms.positions[i]))

        ################ Edit distance in Angstrom here ###################
        if distance_ab > distance_cutoff:
            atoms_to_delete.append(i)

    del atoms[atoms_to_delete]

    return atoms

def transpose(periodic,cluster, start, stop, centre_periodic, centre_cluster, file_name):
    ''' Returns a ChemShell cluster representation of a periodic model 


    Parameters:
    periodic: Atoms object
             Input periodic model of interest
    cluster: Atoms object
             Input cluster model of interest 
    start: Integer 
             Starting index of atoms to copy into cluster
    stop:  Integer
              End index of atoms to copy into cluster
    centre_periodic: Integer
              Centre of periodic model
    centre_cluster: Integer 
              Centre of cluster model
    name: String 
              Name of new model to save 

    '''

    from ase.io import write

    # Start by specifying which atoms in periodic cluster you want to move to a cluster
    to_move = periodic[start:stop]
    #Define central atom in periodic for projection in cluster
    ref_old = periodic[centre_periodic].position
    for i in to_move:
        i.position -= ref_old
    #Define central atom in cluster for reference with periodic
    ref_new = cluster[centre_cluster].position
    for i in to_move:
        i.position += ref_new

    cluster = cluster + to_move
    scaled = str(cluster.get_scaled_positions())
    write = write(file_name, cluster)

def cif2labelpun(frag, charge_dict, in_fname, out_fname, origin_atom, partition_mode='radius',
             radius=5.0, cell_dim=(1, 1, 1)):
    ''' Returns a .pun ChemShell file with region labelled atoms.
        Inherits lattice parameters from Atoms object. Also supports
        radius based partitioning.

        Requires ChemShell.
        May be moved to ChemShell, but currently put here for collaboration work.

    Parameters:
    charge_dict: dict
             Dictionary for charges of each atom
    in_fname: string
             Input file name of periodic structure
    out_fname: string
             Pre-fix filename for .pun output
    origin_atom: integer
              Future functionality for shifting the origin atom
    partition_mode: string
              Accepts 'radius' or 'unit_cell' paritioning modes.
    cell_dim: tuple
              Contains the periodicity of the unit cell QM region
    '''

    import numpy as np
    from chemsh import *

    frag = chemsh.io.tools.convert_atoms_to_frag(read(in_fname), connect_mode='ionic')
    frag.addCharges(charge_dict)

    print(f"Vectors: {frag.cell.a, frag.cell.b, frag.cell.c, frag.cell.alpha, frag.cell.beta, frag.cell.gamma}")
    print(f"coords: {frag.coords}")
    print(f"names: {frag.names}")

    cluster = frag.construct_cluster(radius_cluster=cluster_r, origin_atom=0, adjust_charge=adjust_charge,
                                     radius_active=active_r, bq_margin=bq_margin, bq_density=bq_density,
                                     bq_layer=12.0)

    if partition_mode == 'unit_cell':
        qm_region_coords = expand_cell(frag, cell_dimensions=cell_dim)
        qm_region = match_cell(cluster.coords, qm_region_coords)
    if partition_mode == 'radius':
        qm_region = radius_qm_region(cluster.coords, radius)

    partitioned_cluster = chemsh.cluster.partition.partition(cluster, qm_region=qm_region,
                                                             origin=np.array([0, 0, 0]), cutoff_boundary=4.0,
                                                             interface_exclude=["O"], qmmm_interface='explicit',
                                                             radius_active=20.0)
    # Saving cluster
    partitioned_cluster.save(out_fname + '.pun', 'pun')
    partitioned_cluster.save(out_fname + '.xyz', 'xyz')

    cluster.save('cluster.pun', 'pun')
    frag.save('frag.pun', 'pun')

def expand_cell(frag, cell_dimensions=(1, 1, 1)):
    import numpy as np
    import copy
    from math import ceil
    from ase.units import Bohr

    # Takes the input cell coordinates and returns an expanded cell to
    # designate the QM region. Presumes centering on atom 0.

    vectors = frag.cell.vectors

    print(f"CLEAN COORDS {frag.coords}")

    #    frag.coords[1][0]=(0.5*(vectors[0,0] + vectors[0,1] + vectors[0,2]))
    #    frag.coords[1][1]=(0.5*(vectors[1,0] + vectors[1,1] + vectors[1,2]))
    #    frag.coords[1][2]=(0.5*(vectors[2,0] + vectors[2,1] + vectors[2,2]))

    frag_comp = copy.deepcopy(frag)
    print(f"VECTORS {frag.cell.vectors}")

    if cell_dimensions[0] % 2 == 0:
        imin = -cell_dimensions[0] / 2
        imax = cell_dimensions[0] / 2 + 1
    else:
        imin = -ceil(cell_dimensions[0] / 2)
        imax = ceil(cell_dimensions[0] / 2)

    if cell_dimensions[1] % 2 == 0:
        jmin = -cell_dimensions[1] / 2
        jmax = cell_dimensions[1] + 1
    else:
        jmin = -ceil(cell_dimensions[1] / 2)
        jmax = ceil(cell_dimensions[1] / 2)

    if cell_dimensions[2] % 2 == 0:
        kmin = -cell_dimensions[2] / 2
        kmax = cell_dimensions[2] + 1
    else:
        kmin = -ceil(cell_dimensions[2] / 2)
        kmax = ceil(cell_dimensions[2] / 2)

    for i in np.arange(imin, imax):
        for j in np.arange(jmin, jmax):
            for k in np.arange(kmin, kmax):
                vec = np.matmul(np.array([i, j, k]), vectors)

                frag_comp.coords = np.vstack((frag_comp.coords, frag.coords + vec))

    for i in frag_comp.coords:
        print(i * Bohr)

    frag_comp.save('qmregion.xyz', 'xyz')

    return frag_comp.coords


def match_cell(cluster_coords, frag_coords):
    import numpy as np

    # A slightly convoluted bit of list comprehension outputs a list of overlapping numpy arrays.
    # The [0][0] indexes after the np.where statement are needed to extract the
    # array from the output tuples then the first element from the np.arrays.

    # If statement needed to filter out empty arrays.

    ovlp = [np.where(np.all(np.isclose(cluster_coords, coord, atol=3e-3, rtol=3e-3), axis=1))[0][0] for coord in
            frag_coords if
            np.where(np.all(np.isclose(cluster_coords, coord, atol=3e-3, rtol=3e-3), axis=1))[0].size != 0]

    # Remove duplicates
    ovlp = list(set(ovlp))

    return ovlp


def radius_qm_region(cluster_coords, radius):
    import numpy as np
    # Cuts the cluster using a simple radius model from the origin.
    # Assume origin is 0 for now.

    dist_from_0 = np.linalg.norm(cluster_coords, axis=-1)

    qm_region = [x for x in np.where(dist_from_0 < radius)[0]]

    return qm_region

