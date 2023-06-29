"""

def pre_neb_aims(initial,
                 final,
                 hpc="hawk",
                 basis_set ='light',
                 filename="last_predicted_path.traj"):


    '''
    This function performs a preliminary NEB calculation on user-provided
    structures using ML-NEB. If the calculation does not converge within
    75 steps, it is terminated. Minimum Energy Path energy landscape is
    examined for occurence of multiple local maxima and if detected - geometry
    optimisations on local minima are performed.

    The optimised structures can be used as alternative start/end points for
    further calculations making the NEB calculation easier to converge.

    Parameters:
    hpc: string
        'hawk', 'isambard', 'archer' see carmm.run.aims_path.set_aims_command
    basis_set: string
        'light', 'tight' etc., see carmm.run.aims_path.set_aims_command
    filename: string
        Name of a file containing an unconverged NEB Minimum Energy Path.
        Default is 'last_predicted_path.traj' for CatLearn MLNEB.
    initial: Atoms object
        Starting geometry of a NEB calculation.
    final: Atoms object
        End geometry of a NEB calculation.
    '''
    import os

    if not os.path.exists(filename):
        from ase.io import read
        from catlearn.optimize.mlneb import MLNEB

        # Set the environment parameters
        from carmm.run.aims_path import set_aims_command
        set_aims_command(hpc=hpc, basis_set=basis_set)

        # your settings go here
        def my_calc():
            # New method that gives a default calculator
            from carmm.run.aims_calculator import get_aims_calculator
            return get_aims_calculator(dimensions=2)

        from carmm.build.neb.ilm import neb_identify_local_minima
        from carmm.build.neb.ilm import multiple_local_extrema

        # Desired number of images including start and end point
        # Enough to show energy landscape of Minimum Energy Path
        n = 15

        calculator = my_calc()

        # Setup the Catlearn object for MLNEB
        neb_catlearn = MLNEB(start=initial,
                             end=final,
                             ase_calc=calculator,
                             n_images=n,
                             interpolation='idpp', restart=False)

        # Run the NEB optimisation. Adjust fmax to desired convergence criteria,
        # usually 0.01 ev/A. Max steps set to 75 for preliminary study.
        # MLNEB serial part is quick below 100 structures

        neb_catlearn.run(fmax=0.01,
                         trajectory='ML-NEB.traj',
                         full_output=False,
                         steps=75)

    if multiple_local_extrema(filename=filename) is True:
        print("Multiple extrema detected in the predicted Minimum Energy Path.")
        print("Local minima will be identified and optimised")
        atoms_list, indices = neb_identify_local_minima(filename=filename)

        print(len(atoms_list), "minima detected. Performing geometry optimisations.")

        from ase.optimize import BFGS
        from carmm.run.aims_path import set_aims_command
        from carmm.run.aims_calculator import get_aims_calculator
        set_aims_command(hpc=hpc, basis_set=basis_set)

        x = 0
        for atoms in atoms_list:
            id = indices[x]
            atoms.calc = get_aims_calculator(2, k_grid=(3, 3, 1))
            opt = BFGS(atoms,
                       restart="min_"+str(id)+".pckl",
                       trajectory="min_"+str(id)+".traj")
            opt.run(fmax=0.01)
            x = x+1

        print("Geometry optimisations completed.")
        print("Please consider the structures as alternative start/end points.")

    else:
        print("No multiple extrema detected in the predicted Minimum Energy Path.")

"""