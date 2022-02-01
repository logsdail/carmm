def set_aims_command(hpc='hawk', basis_set='light', defaults=2010):
    '''
    Choose supercomputer and basis_set to obtain FHI-aims run command.
    Can be useful to e.g perform a calculation with a larger basis set
    after a geometry optimisation.

    Parameters:
    hpc: String
        Name of the HPC facility where the jobs is being run
        Options: 'hawk', 'isambard', 'archer', 'young' 
    basis_set: String
        Name of basis set for FHI-aims
        Options: 'light', 'intermediate', 'tight', 'really_tight' etc.
    defaults: int
        Either 2010 or 2020 referring to the default species basis sets
        that come with new FHI-aims release, which adhere to the year 2010
         or 2020 standard. Old 2010 value by default to avoid disruption
         for users.
    '''
    import os



    mpirun = "time mpirun -np $SLURM_NTASKS "
    aprun = "time aprun -n $NPROCS "
    gerun = "gerun "
    srun = "srun --cpu-bind=cores --distribution=block:block --hint=nomultithread "
    executable = "bin/aims.$VERSION.scalapack.mpi.x"
    
    species = "species_defaults/" + "defaults_" + str(defaults) + "/" + basis_set

    if hpc.lower() == 'hawk':
        fhi_aims_directory="/apps/local/projects/scw1057/software/fhi-aims/"
        preamble = mpirun
    elif hpc.lower() == 'isambard':
        fhi_aims_directory="/home/ca-alogsdail/fhi-aims-gnu/"
        preamble = aprun
    #elif hpc.lower() == 'archer': # Retired Jul 2021
    elif hpc.lower() == 'archer2':
        fhi_aims_directory="/work/e05/e05-files-log/shared/software/fhi-aims/"
        preamble = srun
    #elif hpc.lower() == 'thomas': # Retired Oct 2020
    elif hpc.lower() == 'young':
        fhi_aims_directory="/home/mmm0170/Software/fhi-aims/"
        preamble = gerun
    else:
        raise Exception("Inappropriate HPC facility: " + hpc + "is not recognised")

    os.environ["ASE_AIMS_COMMAND"]= preamble + fhi_aims_directory + executable
    os.environ["AIMS_SPECIES_DIR"] = fhi_aims_directory + species

