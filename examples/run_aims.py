#!/usr/bin/env python3

'''
This modules tests aims calculator on different machines
'''

def test_run_aims():
    from carmm.run.aims_path import set_aims_command
    import os

    expected_paths = {
        'hawk': 'time srun --nodes=$SLURM_NNODES --ntasks=$SLURM_NTASKS -d mpirun /apps/local/projects/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'hawk-amd': 'time srun --nodes=$SLURM_NNODES --ntasks=$SLURM_NTASKS -d mpirun /apps/local/projects/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'isambard': 'time aprun -n $NPROCS /home/ca-alogsdail/fhi-aims-gnu/bin/aims.$VERSION.scalapack.mpi.x',
        'archer2': 'srun --cpu-bind=cores --distribution=block:block --hint=nomultithread  /work/e05/e05-files-log/shared/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'young': 'gerun  /home/mmm0170/Software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'aws': 'time srun --mpi=pmi2 --hint=nomultithread --distribution=block:block  apptainer exec /shared/logsdail_group/sing/mkl_aims_2.sif bash /shared/logsdail_group/sing/sing_fhiaims_script.sh $@'
    }

    expected_paths_taskfarm = {
        'hawk': "time srun --nodes=1 --ntasks=40 -d mpirun /apps/local/projects/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x",
        'hawk-amd': 'time srun --nodes=1 --ntasks=64 -d mpirun /apps/local/projects/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'isambard': '',
        'archer2': 'srun --cpu-bind=cores --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=128 /work/e05/e05-files-log/shared/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'young': '',
        'aws': 'time srun --mpi=pmi2 --hint=nomultithread --distribution=block:block --nodes=1 --ntasks=72 apptainer exec /shared/logsdail_group/sing/mkl_aims_2.sif bash /shared/logsdail_group/sing/sing_fhiaims_script.sh $@'
    }

    for hpc in ['hawk', 'hawk-amd', 'isambard', 'archer2', 'young', 'aws']:
        '''Assign the executable command based on HPC'''
        set_aims_command(hpc)

        assert os.environ['ASE_AIMS_COMMAND'] == expected_paths[hpc], f"Path incorrect on {hpc}: {expected_paths[hpc]}\n" \
                                                                      f"Currently: {os.environ['ASE_AIMS_COMMAND']}"

        if hpc in ["hawk", 'hawk-amd', 'archer2', 'aws']:
            set_aims_command(hpc, nodes_per_instance=1)
            assert os.environ['ASE_AIMS_COMMAND'] == expected_paths_taskfarm[hpc], f"Path incorrect on {hpc}: {expected_paths_taskfarm[hpc]}\n" \
                                                                                   f"Currently: {os.environ['ASE_AIMS_COMMAND']}"

    from ase.calculators.aims import Aims
    from carmm.run.aims_calculator import get_aims_and_sockets_calculator
    from carmm.utils.python_env_check import ase_env_check

    for state in range(4):
        #fhi_calc = get_aims_calculator(state)
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(state, verbose=True)

        # Assertion test that the correct calculators are being set
        if ase_env_check('3.22.0'):
            assert(type(sockets_calc.launch_client.calc) == Aims)
        else:
            assert (type(sockets_calc.calc) == Aims)

test_run_aims()
