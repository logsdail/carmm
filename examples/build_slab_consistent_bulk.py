
'''
Testing the build_slab_consistent_bulk functionality
'''

def test_build_slab_consistent_bulk():
    from carmm.build.slab_consistent_bulk_generator import bulk_identifier
    from ase.io import read
    from ase.formula import Formula
    slab = read('data/slab_consistent_bulk/geometry_CoO_111_slab.in')
    new_bulk = bulk_identifier(slab)
    emp_formula = new_bulk.get_chemical_formula(mode='hill', empirical=True)
    total_atoms = 0
    for sym, stoi in Formula(emp_formula).count().items():
        total_atoms += stoi
    formula_units = len(new_bulk)/total_atoms
    assert emp_formula=='CoO'
    assert total_atoms==2
    assert formula_units==3.0


test_build_slab_consistent_bulk()


