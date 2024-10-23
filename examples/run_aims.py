#!/usr/bin/env python3

"""
Tests aims calculator on different machines
"""


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

        assert os.environ['ASE_AIMS_COMMAND'] == expected_paths[
            hpc], f"Path incorrect on {hpc}: {expected_paths[hpc]}\n" \
                  f"Currently: {os.environ['ASE_AIMS_COMMAND']}"

        if hpc in ["hawk", 'hawk-amd', 'archer2', 'aws']:
            set_aims_command(hpc, nodes_per_instance=1)
            assert os.environ['ASE_AIMS_COMMAND'] == expected_paths_taskfarm[
                hpc], f"Path incorrect on {hpc}: {expected_paths_taskfarm[hpc]}\n" \
                      f"Currently: {os.environ['ASE_AIMS_COMMAND']}"

    from ase.calculators.aims import Aims
    from carmm.run.aims_calculator import get_aims_and_sockets_calculator, get_aims_calculator
    from carmm.utils.python_env_check import ase_env_check

    for state in range(4):
        # fhi_calc = get_aims_calculator(state)
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=state, verbose=True)

        default_params = {'relativistic': ('atomic_zora', 'scalar'),
                          'xc': 'pbe',
                          'compute_forces': True,
                          }
        if state == 2:
            default_params['use_dipole_correction'] = 'true'
        if state >= 2:
            default_params['k_grid'] = True

        # Assertion test that the correct calculators and default arguments are being set
        if ase_env_check('3.22.0'):
            assert (type(sockets_calc.launch_client.calc) == Aims)
        else:
            assert (type(sockets_calc.calc) == Aims)
            
        params = getattr(fhi_calc, 'parameters')
        assert params['relativistic'] == ('atomic_zora', 'scalar')
        assert params['xc'] == 'pbe'
        assert params['compute_forces'] is True
        if state == 2:
            assert params['use_dipole_correction'] == 'true'
        if state >= 2:
            assert params['k_grid'] is None

        # libxc test
        sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=state, verbose=True,
                                                                 xc='libxc MGGA_X_MBEEF+GGA_C_PBE_SOL')
        params = getattr(fhi_calc, 'parameters')
        assert params['override_warning_libxc'] == 'true'
        assert params['xc'] == 'libxc MGGA_X_MBEEF+GGA_C_PBE_SOL'

    # Test to make sure that we correctly handle scenario when environment variable isn't
    # set in ASE 3.23. This presents issues downstream, so environment must be set 
    # i.e. executable and species directory.
    from unittest import TestCase
    test_get_aims_exception = TestCase()
    if ase_env_check('3.23.0'):
        with test_get_aims_exception.assertRaises(KeyError):
            del os.environ['ASE_AIMS_COMMAND']
            get_aims_calculator(dimensions=0)


test_run_aims()
