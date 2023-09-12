def set_aims_command(hpc='hawk', basis_set='light', defaults=2010, nodes_per_instance=None):
    """
    Choose supercomputer and basis_set to obtain FHI-aims run command.
    Can be useful to e.g. perform a calculation with a larger basis set
    after a geometry optimisation.

    Parameters:
    hpc: String
        Name of the HPC facility where the jobs are being run
        Options: 'hawk', 'hawk-amd', 'isambard', 'archer2', 'young', 'aws'
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
        "hawk": "time srun",
        "hawk-amd": "time srun",
        "isambard": "time aprun",
        "archer2": "srun --cpu-bind=cores --distribution=block:block --hint=nomultithread",
        "young": "gerun",
        "aws": "time srun --mpi=pmi2 --hint=nomultithread --distribution=block:block"
    }

    assert hpc in preamble, "Inappropriate HPC facility: " + hpc + "is not recognised."

    fhi_aims_directory = {
        "hawk": "/apps/local/projects/scw1057/software/fhi-aims/",
        "hawk-amd": "/apps/local/projects/scw1057/software/fhi-aims/",
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

    no_cpus = {"hawk": 40,
              "hawk-amd": 64,
              "archer2": 128,
              "isambard": 64,
              "young":64,
              "aws": 72}


    if nodes_per_instance:
        #TODO: change name to CPU count and make a separate dictionary for non taskfarmed
        cpu_count_setting = {
            "archer2": f"--nodes={nodes_per_instance} --ntasks={int(no_cpus[hpc] * nodes_per_instance)}",
            "hawk": f"--nodes={nodes_per_instance} --ntasks={int(no_cpus[hpc] * nodes_per_instance)} -d mpirun",
            "hawk-amd": f"--nodes={nodes_per_instance} --ntasks={int(no_cpus[hpc] * nodes_per_instance)} -d mpirun",
            "aws": f"--nodes={nodes_per_instance} --ntasks={int(no_cpus[hpc] * nodes_per_instance)}",
            # TODO: add and test isambard and young task-farmed commands
            "isambard": "",
            "young": "",
        }

        if hpc == "aws":
            assert nodes_per_instance == 1, "FHI-aims does not run on more than one node on AWS at present."

        assert hpc in ["archer2", "hawk", "hawk-amd", "aws"], \
            "Only ARCHER2, Hawk and AWS supported for task-farming at the moment."

    else:
        cpu_count_setting = {
            "hawk": f"--nodes=$SLURM_NNODES --ntasks=$SLURM_NTASKS -d mpirun",
            "hawk-amd": f"--nodes=$SLURM_NNODES --ntasks=$SLURM_NTASKS -d mpirun",
            "isambard": "-n $NPROCS",
            "archer2": "",
            "young": "",
            "aws": ""
        }

    os.environ["ASE_AIMS_COMMAND"] = f"{preamble[hpc]} {cpu_count_setting[hpc]} {executable}"