def set_aims_command(hpc='hawk', basis_set='light', defaults=2010, nodes_per_instance=None):
    """
    Choose supercomputer and basis_set to obtain FHI-aims run command.
    Can be useful to e.g. perform a calculation with a larger basis set
    after a geometry optimisation.

    Parameters:
    hpc: String
        Name of the HPC facility where the jobs are being run
        Options: 'hawk', 'isambard', 'archer2', 'young'
    basis_set: String
        Name of basis set for FHI-aims
        Options: 'light', 'intermediate', 'tight', 'really_tight' etc.
    defaults: int
        Either 2010 or 2020 referring to the default species basis sets
        that come with new FHI-aims release, which adhere to the year 2010
         or 2020 standard. Old 2010 value by default to avoid disruption
         for users.
    """
    import os

    preamble_types = {
        "hawk": "time mpirun -np $SLURM_NTASKS ",
        "isambard": "time aprun -n $NPROCS ",
        "archer2": "srun --cpu-bind=cores --distribution=block:block --hint=nomultithread ",
        "young": "gerun "
    }

    hpc = hpc.lower()
    assert hpc in preamble_types, "Inappropriate HPC facility: " + hpc + "is not recognised."

    fhi_aims_directories = {
        "hawk": "/apps/local/projects/scw1057/software/fhi-aims/",
        "isambard": "/home/ca-alogsdail/fhi-aims-gnu/",
        "archer2": "/work/e05/e05-files-log/shared/software/fhi-aims/",
        "young": "/home/mmm0170/Software/fhi-aims/",
    }

    executable = "bin/aims.$VERSION.scalapack.mpi.x"

    species = "species_defaults/" + "defaults_" + str(defaults) + "/" + basis_set

    preamble = preamble_types[hpc]
    fhi_aims_directory = fhi_aims_directories[hpc]

    if nodes_per_instance:
        assert hpc in ["archer2", "hawk"], "Only ARCHER2 and Hawk supported for task-farming at the moment."

        task_farmed_commands = {
            "archer2": "--nodes=" + str(nodes_per_instance) + " --ntasks=" + str(int(128 * nodes_per_instance)) + " ",
            "hawk": "--nodes=" + str(nodes_per_instance) + " --ntasks=" + str(int(40 * nodes_per_instance)) + " ",
            # TODO: add and test isambard and young task-farmed commands
            "isambard": "",
            "young": "",
        }

        preamble += task_farmed_commands[hpc]

    os.environ["ASE_AIMS_COMMAND"] = preamble + fhi_aims_directory + executable
    os.environ["AIMS_SPECIES_DIR"] = fhi_aims_directory + species
