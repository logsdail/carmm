# Author: Igor Kowalec
'''
This example shows how to use the work flow functionality with CARMM to optimise
reactants/products and find transition states in an automated manner using the MACE-MP force field.
'''


def test_run_workflows_ReactMACE():
    '''

    Args:
        model_size: str
            small medium large - size of the mace-mp model to be used

    Returns:

    '''
    from carmm.run.workflows.react_mace import ReactMACE
    from ase.build import molecule, bulk
    from ase.build import surface, add_adsorbate
    from carmm.analyse.forces import is_converged
    from carmm.build.neb.symmetry import sort_z
    from ase.constraints import FixAtoms
    import os
    import math

    '''Work in a dedicated folder'''
    parent_dir = os.getcwd()
    os.makedirs('data/react_mace', exist_ok=True)
    os.chdir('data/react_mace')

    '''Determine calculation input settings'''
    model_size = 'small'
    params = {'model': model_size,
              'dispersion': False,
              'default_dtype': 'float64',
              'device': 'cpu'}

    '''Build the React_Aims object using the set of settings'''
    reactor = ReactMACE(params)

    '''Prepare the Atoms object geometry'''
    '''Calculate adsorbate'''
    adsorbate = molecule('H')

    reactor.filename = "H2"
    ref_adsorbate = molecule("H2")
    ref_adsorbate = reactor.mace_optimise(ref_adsorbate, fmax=0.01, restart=False)
    e_adsorbate = ref_adsorbate.get_potential_energy()

    assert math.sqrt((e_adsorbate - (-6.557687))**2) < 1e-6, "H2 energy with a small dataset should produce -6.557687"

    '''Calculate bulk energy'''
    reactor.filename = f"Pd_bulk_{model_size}"
    bulk_geometry = bulk("Pd", a=3.96)
    bulk_geometry = reactor.mace_optimise(bulk_geometry, fmax=0.01, relax_unit_cell=True, restart=True)
    e_bulk = bulk_geometry.get_potential_energy()
    

    '''Calculate pristine surface slab energy'''
    reactor.filename = f"Pd_111_{model_size}"
    slab = surface(bulk_geometry, (1,1,1), 7, vacuum=20)
    slab = slab.repeat((3,3,1))
    slab = sort_z(slab, diff=1)
    constraint = FixAtoms(indices=[atom.index for atom in slab if atom.tag > 4])
    slab.set_constraint([constraint])
    slab = reactor.mace_optimise(slab, fmax=0.01, restart=True)
    e_pristine = slab.get_potential_energy()
    reactor.verbose = False

    '''Calculate slab with adsorbate'''
    adsorption_sites = [
        {"site": "fcc", "position": 25},
        {"site": "hcp", "position": 26},
        # {"site": "atop", "position": 56}
        ]

    neb_input = []
    for configuration in adsorption_sites:
        reactor.filename = f"Pd_111_H_{configuration['site']}_{model_size}"
        slab_copy = slab.copy()
        add_adsorbate(slab_copy, adsorbate, height=1, position=(slab[configuration['position']].position[:2]))
        slab_with_adsorbate = reactor.mace_optimise(slab_copy, fmax=0.01, restart=True)
        e_slab = slab_with_adsorbate.get_potential_energy()

        neb_input.append(slab_with_adsorbate)

        '''Derive the adsorption energy'''
        e_ads = e_slab - (e_pristine + 0.5 * e_adsorbate)
        print(f"Using {model_size} model, the E_ads at {configuration['site']} site is {e_ads} eV")

    ts = reactor.search_ts_neb(initial=neb_input[0],
                               final=neb_input[1],
                               fmax=0.05,
                               k=0.05,
                               n=7,
                               interpolation="idpp",
                               input_check=0.05,
                               restart=True)

    '''Call relevant calculations'''
    '''Dissociate hydrogen using MACE forcefield'''
    reactor = ReactMACE(params)
    H2 = ref_adsorbate.copy()
    print(H2)
    H_H = ref_adsorbate.copy()
    H_H[0].x += 7
    ts = reactor.search_ts_neb(initial=H2, final=H_H, fmax=0.05, k=0.05, n=5, interpolation="idpp", input_check=0.05,
                               restart=True)

    '''Test reusing manually provided interpolation'''
    ts = reactor.search_ts_neb(initial=H2, final=H_H, fmax=0.05, k=0.05, n=5, interpolation=reactor.interpolation,
                               input_check=0.05, restart=False)

    assert is_converged(reactor.model_optimised, 0.05), \
        '''The structure saved in React_Aims is not converged'''

    os.chdir(parent_dir)

test_run_workflows_ReactMACE()