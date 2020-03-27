def get_aims_command(supercomputer='hawk', basis_set='light'):

    # Choose supercomputer and basis_set to obtain FHI-aims run command.
    # Can be useful to e.g perform a calcaluation with a larger basis set
    # after a geometry optimisation.
    # supercomputers: 'hawk', 'isambard', 'archer'  ## needs Thomas
    # basis_set : 'light', 'tight', 'really_tight' etc.

    if supercomputer =='hawk':
        os.environ["ASE_AIMS_COMMAND"]="mpirun -np $SLURM_NTASKS /home/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x"
        os.environ['AIMS_SPECIES_DIR']="/home/scw1057/software/fhi-aims/species_defaults/"+basis_set

    elif supercomputer =='isambard':
        import os
        command="aprun -n "+os.environ["NPROCS"]+" /home/ca-alogsdail/fhi-aims-gnu/bin/aims."+os.environ["VERSION"]+".scalapack.mpi.x"
        os.environ["ASE_AIMS_COMMAND"]=command
        os.environ["AIMS_SPECIES_DIR"]="/home/ca-alogsdail/fhi-aims-gnu/species_defaults/"+basis_set

    elif supercomputer =='archer':
        command="aprun -n "+os.environ["NPROCS"]+" /home3/e05/e05/ajl340/fhi-aims-src-intel/bin/aims."+os.environ["VERSION"]+".scalapack.mpi.x"
        os.environ["ASE_AIMS_COMMAND"]=command
        os.environ["AIMS_SPECIES_DIR"]="/home3/e05/e05/ajl340/fhi-aims-src-intel/species_defaults/"+basis_set
