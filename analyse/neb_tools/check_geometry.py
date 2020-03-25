from ase import Atoms
from ase.io import read,write
from ase.io.trajectory import Trajectory

def switch_indices(model,A,B):
    if isinstance(model, str) is True:
       model = read(model)

    e=model.get_potential_energy()
    f=model.get_forces()
    new_model=Atoms()

    # index manipulation
    if not isinstance(A, int) and isinstance(B, int):
        raise ValueError
        print("Indices must be integers.")

    else:
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
    t1=Trajectory('switched_i_'+str(A)+'_'+str(B)+'.traj', 'w')
    t1.write(new_model, energy=e, forces=f)
