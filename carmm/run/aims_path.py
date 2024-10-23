def set_aims_command(hpc='hawk', basis_set='light', defaults=2010, nodes_per_instance=None):
    """
    Choose supercomputer and basis_set to obtain FHI-aims run command.
    Can be useful to e.g. perform a calculation with a larger basis set
    after a geometry optimisation.

    Parameters:
    hpc: String
        Name of the HPC facility where the jobs are being run
        Options: 'hawk', 'hawk-amd', 'isambard', 'archer2', 'young', 'aws', 'custom'
        NOTE 1: 'custom' requires the environmental variable "CARMM_AIMS_ROOT_DIRECTORY"
        before running to allow logic of basis set selection, while maintaining
        free choice of basis set folders.
        NOTE 2: 'custom' requires environmental variable ASE_AIMS_COMMAND be set, with
        the desired number of mpi processes and WITHOUT a default output file
    basis_set: String
        Name of basis set for FHI-aims
        Options: 'light', 'intermediate', 'tight', 'really_tight' etc.
    defaults: int
        Either 2010 or 2020 referring to the default species basis sets
        that come with new FHI-aims release, which adhere to the year 2010
         or 2020 standard. Old 2010 value by default to avoid disruption
         for users.
    nodes_per_instance: int, optional
        Number of nodes per separate instance of FHI-aims, when running task-farmed
    """
    import os

    hpc = hpc.lower()

    if hpc == "custom":
        assert "CARMM_AIMS_ROOT_DIRECTORY" in os.environ, \
            "hpc is 'custom' but environmental variable CARMM_AIMS_ROOT_DIRECTORY not specified."

        custom_root_dir = os.environ["CARMM_AIMS_ROOT_DIRECTORY"]
    else:
        custom_root_dir = None


    species = "species_defaults/" + "defaults_" + str(defaults) + "/" + basis_set

    preamble = {
        "hawk": "time srun",
        "hawk-amd": "time srun",
        "isambard": "time aprun",
        "archer2": "srun --cpu-bind=cores --distribution=block:block --hint=nomultithread",
        "young": "gerun",
        "aws": "time srun --mpi=pmi2 --hint=nomultithread --distribution=block:block",
        "custom": ""
    }

    assert hpc in preamble, "Inappropriate HPC facility: " + hpc + "is not recognised."

    fhi_aims_directory = {
        "hawk": "/apps/local/projects/scw1057/software/fhi-aims/",
        "hawk-amd": "/apps/local/projects/scw1057/software/fhi-aims/",
        "isambard": "/home/ca-alogsdail/fhi-aims-gnu/",
        "archer2": "/work/e05/e05-files-log/shared/software/fhi-aims/",
        "young": "/home/mmm0170/Software/fhi-aims/",
        "aws": "/shared/logsdail_group/sing/",
        "custom": custom_root_dir
    }

    executable_d = {"compiled": "bin/aims.$VERSION.scalapack.mpi.x",
                    "apptainer": "apptainer exec " + fhi_aims_directory["aws"] + "mkl_aims_2.sif bash " + \
                                 fhi_aims_directory["aws"] + "sing_fhiaims_script.sh $@"
                    }

    '''Handle compiled and containerized FHIaims versions'''
    if hpc == "aws":
        executable = executable_d["apptainer"]
    elif hpc != "custom":
        executable = fhi_aims_directory[hpc] + executable_d["compiled"]

    """Set the relevant environment variables based on HPC"""
    os.environ["AIMS_SPECIES_DIR"] = fhi_aims_directory[hpc] + species

    # Define the executable command
    if nodes_per_instance:
        # Check validity of task-farming setup before proceeding.
        # Todo: Add Isambard/Young as needed
        assert hpc in ["archer2", "hawk", "hawk-amd", "aws"], \
            "Only ARCHER2, Hawk and AWS supported for task-farming at the moment."

        if hpc == "aws":
            assert nodes_per_instance == 1, "FHI-aims does not run on more than one node on AWS at present."

    if hpc == 'custom':
        assert "ASE_AIMS_COMMAND" in os.environ, \
            "set_aims_command: option hpc is 'custom', but ASE_AIMS_COMMAND not set."
    else:
        # This has a helper function as we need to take different actions
        # if running single or task-farmed calculations
        cpu_command = _get_cpu_command(hpc, nodes_per_instance)

        os.environ["ASE_AIMS_COMMAND"] = f"{preamble[hpc]} {cpu_command} {executable}"

def _get_cpu_command(hpc, nodes_per_instance=None):
    """
    Helper function to return appropriate syntax for cpu settings on each HPC

    Parameters:
        As for set_aims_command
    """

    # This dictionary contains settings related to each HPC infrastructure
    hpc_settings = {
        "hawk": { "cpus_per_node": 40, "cpu_command": f"--nodes=$SLURM_NNODES --ntasks=$SLURM_NTASKS -d mpirun", },
        "hawk-amd": { "cpus_per_node": 64, "cpu_command": f"--nodes=$SLURM_NNODES --ntasks=$SLURM_NTASKS -d mpirun", },
        "isambard": { "cpus_per_node": 64, "cpu_command": f"-n $NPROCS", },
        "young": { "cpus_per_node": 64, "cpu_command": "", },
        "archer2": { "cpus_per_node": 128, "cpu_command": "", },
        "aws": { "cpus_per_node": 72, "cpu_command": "",  }
    }
    
    # This content adds capabilities relating to task-farming.
    # Todo: Extend for Isambard/Young
    if nodes_per_instance:
        hpc_settings["hawk"]["cpu_command_task_farming"] = f"--nodes={nodes_per_instance} --ntasks={int(hpc_settings['hawk']['cpus_per_node'] * nodes_per_instance)} -d mpirun"
        hpc_settings["hawk-amd"]["cpu_command_task_farming"] = f"--nodes={nodes_per_instance} --ntasks={int(hpc_settings['hawk-amd']['cpus_per_node'] * nodes_per_instance)} -d mpirun"
        hpc_settings["archer2"]["cpu_command_task_farming"] = f"--nodes={nodes_per_instance} --ntasks={int(hpc_settings['archer2']['cpus_per_node'] * nodes_per_instance)}"
        hpc_settings["aws"]["cpu_command_task_farming"] = f"--nodes={nodes_per_instance} --ntasks={int(hpc_settings['aws']['cpus_per_node'] * nodes_per_instance)}"

    # Check calculation effiency
    if hpc in ["hawk","hawk-amd"]:
        # Necessary import
        import os
        
        # Placed inside try/except to work when environment variables aren't defined
        try: requested_tasks = int(os.environ["SLURM_NTASKS"])
        except: requested_tasks = hpc_settings[hpc]["cpus_per_node"]

        try: requested_nodes = int(os.environ["SLURM_NNODES"])
        except: requested_nodes = 1

        # Check if using a full node, for efficiency
        if (requested_tasks/requested_nodes) % hpc_settings[hpc]["cpus_per_node"] != 0:
            print("WARNING: You are not using all the CPUs on the requested nodes.")
            print("         Check if you are accidentally underpopulating the nodes.")  

    if nodes_per_instance:
        return hpc_settings[hpc]["cpu_command_task_farming"]
    else:
        return hpc_settings[hpc]["cpu_command"]


    
    
