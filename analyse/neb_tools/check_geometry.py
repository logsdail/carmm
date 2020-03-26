from ase import Atoms
from ase.io import read,write
from ase.io.trajectory import Trajectory

def switch_indices(model, A, B, output_file):
    ## Takes 4 arguments: filename/Atoms object,
    ## indices A,B to switch around
    ## and output_file 'name.traj'

    if isinstance(model, str) is True:
       model = read(model)

    e=model.get_potential_energy()
    f=model.get_forces()
    new_model=Atoms()

    # index manipulation
    if not isinstance(A, int) and isinstance(B, int):
        raise ValueError
        print("Indices must be integers.")

    elif A ==B:
        print("Chosen idices are identical.")

    elif B < A:
            A,B = B,A # make sure indices are in the right order

    for i in range(0,A):
        new_model += model[i]
    new_model +=model[B]
    for i in range(A+1,B):
        new_model += model[i]
    new_model += model[A]
    for i in range(B+1,len(model.get_tags())):
        new_model += model[i]
    # Copy parameters into the new model
    new_model.set_cell(model.get_cell())
    new_model.set_pbc(model.get_pbc())
    new_model.set_constraint(model._get_constraints())

    # writer section
    #write('Model_corrected.traj', new_model)
    t1=Trajectory(str(output_file), 'w')
    t1.write(new_model, energy=e, forces=f)
    t1.close()

def check_interpolation(initial, final, n_max):
    ## interpolates the provided geometries with n_max total images
    ## and checks whether any bond lengths below 0.74 Angstrom exist
    ## saves the intepolation in interpolation.traj

    from ase import io
    import copy
    from ase.neb import NEB
    from software.analyse.Interatomic_distances.analyse_bonds import search_abnormal_bonds

    # pre-requirements
    if isinstance(initial, str) is True:
       initial = read(initial)
    if isinstance(final, str) is True:
       final = read(final)
    if not isinstance(n_max, int):
        raise ValueError
        print("Max number of images must be an integer.")

    # Make a band consisting of 10 images:
    images = [initial]
    images += [initial.copy() for i in range(n_max-2)]
    images += [final]
    neb = NEB(images, climb=True)
    # Interpolate linearly the potisions of the middle images:
    neb.interpolate()

    t2=Trajectory('interpolation.traj', 'w')

    for i in range(0,n_max):
        print("Assessing image",str(i+1)+'.')
        search_abnormal_bonds(images[i])
        t2.write(images[i])

    t2.close()
