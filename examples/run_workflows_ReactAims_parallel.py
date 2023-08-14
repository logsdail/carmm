'''
Author: Igor Kowalec

This example shows how to use the workflow functionality with CARMM to optimise
structures and perform efficient task-farming with concurrent.futures module in an ASE/FHI-aims setup.
Subsequent jobs will be submitted as soon as CPUs become available resulting in time savings
compared to taskfarming in batches.
The job must still run until completion and some resource waste towards the end
of the calculations is expected.
'''


def test_run_workflows_ReactAims_parallel():
    from carmm.run.workflows.react import ReactAims
    import csv
    import concurrent.futures
    from ase.build import molecule
    import os

    '''Work in a dedicated folder'''
    parent_dir = os.getcwd()
    os.makedirs("data/react", exist_ok=True)
    os.chdir("data/react")

    '''Run the ASE/FHI-aims/CatLearn transition state search on the hpc facility'''
    TOTAL_NODES = 4
    NODES_PER_INSTANCE = 1
    PARAMS = {"xc":'pbe',
             "spin":'none',
            }
    BASIS_SET = "light"
    HPC = "archer2"

    '''Input previously generated and calculated ASE Atoms object - pristine surface'''
    molecules = [f"C{n}H{2*n+2}" for n in range(2,4)]*10
    sampled_structures = [[f"Test_{i}", molecule(molecules[i])] for i in range(20)]

    '''Prepare the tasks in the form of a generator'''
    def condensed_generator(n):
        params = PARAMS.copy()
        reactor= ReactAims(params=params,
                           basis_set=BASIS_SET,
                           hpc=HPC,
                           filename=sampled_structures[n][0],
                           nodes_per_instance=NODES_PER_INSTANCE)

        reactor.aims_optimise(sampled_structures[n][1],
                                              fmax=0.01,
                                              restart=True)
        e = reactor.model_optimised.get_potential_energy()
        configuration_info = [reactor.filename, e]

        return configuration_info


    # List of task parameters
    # More parameters can be provided in a tuple
    task_parameters_list = [(n) for n in range(len(sampled_structures))]

    # Number of tasks to execute concurrently
    concurrent_tasks = int(TOTAL_NODES / NODES_PER_INSTANCE)

    # Using ThreadPoolExecutor for concurrent execution (you can also use ProcessPoolExecutor)
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrent_tasks) as executor:
        with open("config_info.csv", "w", newline="") as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=",")
            # Submit tasks to the executor and collect the futures
            futures = [executor.submit(condensed_generator, n) for n in task_parameters_list]
            # You can retrieve results as they complete
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                # Process the result if needed
                # Write the entire result list to the CSV
                csvwriter.writerow(result)

    '''Return to parent directory'''
    os.chdir(parent_dir)

test_run_workflows_ReactAims_parallel()