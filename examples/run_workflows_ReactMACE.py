'''
This example shows how to use the work flow functionality with CARMM to optimise
reactants/products and find transition states in an automated manner using the MACE-MP force field.
'''


def test_run_workflows_ReactAims():
    from carmm.run.workflows.react_mace import ReactMACE
    from ase.build import molecule, bulk
    from carmm.analyse.forces import is_converged
    import os

    '''Work in a dedicated folder'''
    parent_dir = os.getcwd()
    os.makedirs('data/react_mace', exist_ok=True)
    os.chdir('data/react_mace')

    '''Determine calculation input settings'''
    params = {'model':'small',
              'default_dtype':'float64',
              'device':'cpu'}

    '''Prepare the Atoms object geometry'''
    atoms = molecule('H2')

    '''Build the React_Aims object using the set of settings'''
    reactor = ReactMACE(params, dry_run=True)

    '''Call relevant calculations'''
    '''The below has been previously calculated and data is retrieved from saved trajectories'''
    model_optimised = reactor.mace_optimise(atoms, fmax=0.05, restart=True)

    assert is_converged(reactor.model_optimised, 0.01), \
        '''The structure saved in React_Aims is not converged'''


test_run_workflows_ReactAims()