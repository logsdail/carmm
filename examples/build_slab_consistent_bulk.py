
'''
Testing the build_slab_consistent_bulk functionality
'''

def test_build_slab_consistent_bulk():
    # importing libraries inside the function
    from carmm.build.slab_consistent_bulk_generator import bulk_identifier
    from ase.formula import Formula
    from ase.atoms import Atoms
    from ase.build import bulk
    import numpy
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core.surface import SlabGenerator
    crys = bulk('MgO', crystalstructure='rocksalt', cubic=True, a=4.21)
    structure = AseAtomsAdaptor.get_structure(crys)
    slabgen = SlabGenerator(structure, miller_index=(1, 1, 1),
                            min_slab_size=10,
                            min_vacuum_size=20,
                            center_slab=True, in_unit_planes=True, lll_reduce=True)
    slabs = slabgen.get_slabs(ftol=0.001, symmetrize=False)
    slab = AseAtomsAdaptor.get_atoms(slabs[0].get_orthogonal_c_slab())
    slab = Atoms(slab)
    new_bulk = bulk_identifier(slab)
    emp_formula = new_bulk.get_chemical_formula(mode='hill', empirical=True)
    total_atoms = 0

    for sym, stoi in Formula(emp_formula).count().items():
        total_atoms += stoi
    formula_units = len(new_bulk)/total_atoms

    assert emp_formula=='MgO'
    assert total_atoms==2
    assert formula_units==3.0


test_build_slab_consistent_bulk()


