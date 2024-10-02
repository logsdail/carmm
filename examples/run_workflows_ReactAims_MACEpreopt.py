'''
This example shows how to use the work flow functionality with CARMM to optimise
reactants/products and find transition states in an automated manner using FHI-aims.
@Igor: Can you update with any further information?
'''

def test_run_workflows_ReactAims_MACE_preopt():
    from carmm.run.workflows.react import ReactAims
    from carmm.run.workflows.react_mace import ReactMACE
    from ase.build import molecule, bulk
    from carmm.analyse.forces import is_converged
    import os
    import numpy as np
    from ase import Atoms
    from ase.optimize import FIRE
    from ase.build import surface, add_adsorbate
    from ase.constraints import FixAtoms

    '''Work in a dedicated folder'''
    parent_dir = os.getcwd()
    os.makedirs("data/react_preopt", exist_ok=True)
    os.chdir("data/react_preopt")

    '''Determine calculation input settings'''
    params = {"xc":"pbe"}
    basis_set = "light"
    hpc = "hawk"

    mace_params = {'model': 'small',
                   'dispersion': False,
                   'default_dtype': 'float64',
                   'device': 'cpu'}
    
    '''Prepare the Atoms object geometry'''
    atoms = molecule("H2")

    '''Build the React_Aims object using the set of settings'''
    reactor = ReactAims(params, basis_set, hpc, dry_run=True)
    reactor.MaceReact_Preoptimiser = ReactMACE(mace_params)

    '''Call relevant calculations'''
    '''The below has been previously calculated and data is retrieved from saved trajectories'''
    model_optimised, _ = reactor.aims_optimise(atoms, fmax=0.01, restart=True, optimiser=FIRE, mace_preopt=True)
    '''I do not like this test - it tests EMT more than it does the workflow'''
    assert is_converged(reactor.model_optimised, 0.01), \
        "Model is not optimised after preoptimisation"

    '''Create a reaction pathway'''
    atoms1 = atoms.copy()
    atoms1[1].x += 8
    atoms2 = atoms.copy()
    '''Tests each flavour of ML-NEB'''
    transition_state = reactor.search_ts(atoms1, atoms2, 0.05, 0.03, n=6, input_check=0.01, mace_preopt="tspath", restart=False)
    print(len(reactor.mace_interpolation))
    assert is_converged(reactor.mace_interpolation[3], 0.03), \
        "MACE pre-optimised NEB flavour 2 is not converged"
    transition_state = reactor.search_ts(atoms1, atoms2, 0.05, 0.03, n=6, input_check=0.01, mace_preopt="fullpath", restart=False)
    assert is_converged(reactor.mace_interpolation[3], 0.03), \
        "MACE pre-optimised NEB flavour 1 is not converged"
    transition_state = reactor.search_ts(atoms1, atoms2, 0.05, 0.03, n=6, input_check=0.01, mace_preopt=None, restart=False)
    assert isinstance(reactor.interpolation, str), \
        "ML-NEB is somehow running even though not requested"

    '''Return to parent directory'''
    os.chdir(parent_dir)


test_run_workflows_ReactAims_MACE_preopt()

