def set_aims_command(supercomputer='hawk', basis_set='light'):
    import os
    # Choose supercomputer and basis_set to obtain FHI-aims run command.
    # Can be useful to e.g perform a calcaluation with a larger basis set
    # after a geometry optimisation.
    # supercomputers: 'hawk', 'isambard', 'archer'  ## needs Thomas
    # basis_set : 'light', 'tight', 'really_tight' etc.

    if supercomputer =='hawk':
        os.environ["ASE_AIMS_COMMAND"]="time mpirun -np $SLURM_NTASKS /home/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x"
        os.environ['AIMS_SPECIES_DIR']="/home/scw1057/software/fhi-aims/species_defaults/"+basis_set

    elif supercomputer =='isambard':
        os.environ["ASE_AIMS_COMMAND"]="time aprun -n $NPROCS /home/ca-alogsdail/fhi-aims-gnu/bin/aims.$VERSION.scalapack.mpi.x"
        os.environ["AIMS_SPECIES_DIR"]="/home/ca-alogsdail/fhi-aims-gnu/species_defaults/"+basis_set

    elif supercomputer =='archer':
        os.environ["ASE_AIMS_COMMAND"]="time aprun -n $NPROCS /home3/e05/e05/ajl340/fhi-aims-src-intel/bin/aims.$VERSION.scalapack.mpi.x"
        os.environ["AIMS_SPECIES_DIR"]="/home3/e05/e05/ajl340/fhi-aims-src-intel/species_defaults/"+basis_set
