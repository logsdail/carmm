def set_aims_command(hpc='hawk', basis_set='light', defaults=2010, nodes_per_instance=None):
    """
    Choose supercomputer and basis_set to obtain FHI-aims run command.
    Can be useful to e.g. perform a calculation with a larger basis set
    after a geometry optimisation.

    Parameters:
    hpc: String
        Name of the HPC facility where the jobs are being run
        Options: 'hawk', 'isambard', 'archer2', 'young', 'aws'
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

    hpc = hpc.lower()

    species = "species_defaults/" + "defaults_" + str(defaults) + "/" + basis_set

    preamble = {
        "hawk": "time mpirun -np $SLURM_NTASKS ",
        "isambard": "time aprun -n $NPROCS ",
        "archer2": "srun --cpu-bind=cores --distribution=block:block --hint=nomultithread ",
        "young": "gerun ",
        "aws": "time srun --mpi=pmi2 --hint=nomultithread --distribution=block:block "
    }

    assert hpc in preamble, "Inappropriate HPC facility: " + hpc + "is not recognised."

    fhi_aims_directory = {
        "hawk": "/apps/local/projects/scw1057/software/fhi-aims/",
        "isambard": "/home/ca-alogsdail/fhi-aims-gnu/",
        "archer2": "/work/e05/e05-files-log/shared/software/fhi-aims/",
        "young": "/home/mmm0170/Software/fhi-aims/",
        "aws": "/shared/logsdail_group/sing/",
    }

    executable_d = {"compiled": "bin/aims.$VERSION.scalapack.mpi.x",
                    "apptainer": "apptainer exec " + fhi_aims_directory["aws"] + "mkl_aims_2.sif bash " + \
                                 fhi_aims_directory["aws"] + "sing_fhiaims_script.sh $@"
                    }

    '''Handle compiled and containerized FHIaims versions'''
    if hpc == "aws":
        executable = executable_d["apptainer"]
    else:
        executable = fhi_aims_directory[hpc] + executable_d["compiled"]

    """Set the relevant environment variables based on HPC"""
    os.environ["AIMS_SPECIES_DIR"] = fhi_aims_directory[hpc] + species
    if nodes_per_instance:
        task_farmed_commands = {
            "archer2": "--nodes=" + str(nodes_per_instance) + " --ntasks=" + str(int(128 * nodes_per_instance)) + " ",
            "hawk": "--nodes=" + str(nodes_per_instance) + " --ntasks=" + str(int(40 * nodes_per_instance)) + " ",
            # TODO: add and test isambard and young task-farmed commands
            "isambard": "",
            "young": "",
            "aws": "--nodes=" + str(nodes_per_instance) + " --ntasks=" + str(int(72 * nodes_per_instance)) +" ",
        }

        assert hpc in ["archer2", "hawk", "aws"], "Only ARCHER2, Hawk and AWS supported for task-farming at the moment."
        os.environ["ASE_AIMS_COMMAND"] = preamble[hpc] + task_farmed_commands[hpc] + executable
    else:
        os.environ["ASE_AIMS_COMMAND"] = preamble[hpc] + executable
