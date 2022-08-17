#!/usr/bin/env python3

'''
This modules tests aims calculator on different machines - Hawk, Isambard and archer.
The dimensions from 0 to 3 test the different aims setting as described in
carmm.run.get_aims_calculator
TODO: rework the assertion to actually test the aims_path output - include expected string values
'''

def test_run_aims():
    from carmm.run.aims_path import set_aims_command
    import os

    expected_paths = {
        'hawk': 'time mpirun -np $SLURM_NTASKS /apps/local/projects/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'hawk-amd': 'time mpirun -np $SLURM_NTASKS /apps/local/projects/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'isambard': 'time aprun -n $NPROCS /home/ca-alogsdail/fhi-aims-gnu/bin/aims.$VERSION.scalapack.mpi.x',
        'archer2': 'srun --cpu-bind=cores --distribution=block:block --hint=nomultithread /work/e05/e05-files-log/shared/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'young': 'gerun /home/mmm0170/Software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'aws': 'time srun --mpi=pmi2 --hint=nomultithread --distribution=block:block apptainer exec /shared/logsdail_group/sing/mkl_aims_2.sif bash /shared/logsdail_group/sing/sing_fhiaims_script.sh $@'
    }

    expected_paths_taskfarm = {
        'hawk': "time mpirun -np $SLURM_NTASKS --nodes=1 --ntasks=40 /apps/local/projects/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x",
        'hawk-amd': 'time mpirun -np $SLURM_NTASKS --nodes=1 --ntasks=64 /apps/local/projects/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'isambard': '',
        'archer2': 'srun --cpu-bind=cores --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=128 /work/e05/e05-files-log/shared/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x',
        'young': '',
        'aws': 'time srun --mpi=pmi2 --hint=nomultithread --distribution=block:block --nodes=1 --ntasks=72 apptainer exec /shared/logsdail_group/sing/mkl_aims_2.sif bash /shared/logsdail_group/sing/sing_fhiaims_script.sh $@'
    }

    for hpc in ['hawk', 'hawk-amd', 'isambard', 'archer2', 'young', 'aws']:
        set_aims_command(hpc)
        assert os.environ['ASE_AIMS_COMMAND'] == expected_paths[hpc]

        if hpc in ["hawk", 'hawk-amd', 'archer2', 'aws']:
            set_aims_command(hpc, nodes_per_instance=1)
            assert os.environ['ASE_AIMS_COMMAND'] == expected_paths_taskfarm[hpc]

    import ase # Necessary to check the code version, as socket functionality has changed
    ase_major_version = int(ase.__version__.split(".")[0])
    ase_minor_version = int(ase.__version__.split(".")[1])

    from ase.calculators.aims import Aims
    from carmm.run.aims_calculator import get_aims_and_sockets_calculator

    for state in range(4):
        #fhi_calc = get_aims_calculator(state)
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(state, verbose=True)

        # Assertion test that the correct calculators are being set
        # ASE version 3.21 or earlier
        if ase_major_version <= 3 and ase_minor_version <= 21:
            assert (type(sockets_calc.calc) == Aims)
        else:
        # ASE Version 3.22 or later
            assert(type(sockets_calc.launch_client.calc) == Aims)

test_run_aims()
