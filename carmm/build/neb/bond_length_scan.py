'''
Created on Fri 19/06/2020

@author: Igor Kowalec, David Willock
'''

def dissociation(atoms, i1, i2, step_size=0.05, n_steps=20, traj_prefix="mep", fmax=0.05):
    '''
    This function is a tool for investigating atom dissociation from a molecule.
    Bond length of interest is fixed and is increased by step_size in each iteration.
    This aims to help with obtaining activation energies for surface calculations, e.g. hydrogenation,
    where metastability of optimal starting position is often low and thus hard to obtain.

    Args:
        atoms: Atoms object
        i1: int
            Index of atom remaining as part of a molecule
        i2: int
            Index of atom dissociating from molecule
        step_size: float
            Distance moved during dissociation in Angstrom per iteration
        n_steps: int
            Number of steps to calculate
        traj_prefix: str
            Naming scheme that will be followed for all geometry optimisation .traj files
        fmax: float
            Force convergence criterion for geometry optimisation, default is 0.01 eV/Angstrom
    '''

    from ase.io.trajectory import Trajectory
    from ase.optimize import BFGS
    from ase.constraints import FixBondLength
    import copy

    # Helper function Atom-Atom distance
    def atom_dist(atoms, i1, i2):
        import numpy as np
        pos_diff = atoms[i1].position - atoms[i2].position
        distance = np.linalg.norm(pos_diff)

        return distance


    # retrieve initial atom - atom distance
    initial_dist = atom_dist(atoms, i1, i2)
    if not atoms.constraints:
        initial_constraint = None
    else:
        initial_constraint = copy.deepcopy(atoms.constraints)

    traj_counter = 0

    # Write out the initial following naming scheme
    traj_initial = Trajectory(str(traj_prefix) + str(traj_counter) + ".traj", 'w')
    traj_initial.write(atoms)

    for i in range(1, n_steps+1):
        traj_counter += 1

        #remove previous contraints and set up new ones
        atoms.set_constraint()
        new_constraint = FixBondLength(i1, i2)
        atoms.set_distance(i1, i2, (initial_dist + i * step_size), fix=0)

        if initial_constraint is not None:
            atoms.set_constraint([initial_constraint, new_constraint])
        else:
            atoms.set_constraint(new_constraint)
        opt = BFGS(atoms, trajectory=traj_prefix+str(traj_counter)+".traj",
                   restart=traj_prefix+str(traj_counter)+".pckl")
        opt.run(fmax=fmax)

    print("Have a nice day.")

def test_dissociation():
    from ase.build import fcc100, add_adsorbate
    from ase.calculators.emt import EMT
    from ase import Atoms
    slab = fcc100("Cu", (3,3,4), vacuum=10, periodic=True)
    Au_adatoms = Atoms("2Au", positions=[(0,0,0), (2,0,0)])

    add_adsorbate(slab, Au_adatoms, 2, position=(slab[18].x, slab[18].y))
    slab.set_calculator(EMT())

    # start with an optimised structure
    from ase.optimize import BFGS
    opt = BFGS(slab)
    opt.run(fmax=0.05)

    # move gold atoms away from each other
    dissociation(slab, 36, 37, step_size=0.1)

test_dissociation()