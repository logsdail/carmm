#!/usr/bin/env python3
'''
This example shows usage of symmetry operation analysis for FCC lattice.
So far, functionality is implemented only for low index surfaces:
(111), (110) and (100)

These come in handy when setting up NEBs.
'''

def test_build_symmetry_mirror():

    ###### EXAMPLE OF USE - MIRROR AND TRANSLATION ########
    #from ase.visualize import view
    from data.model_gen import get_example_slab as slab
    from carmm.build.neb.symmetry import mirror, translation, rotate_fcc

    for surface_facet in ['100', '110']:
        # Toy model of CO2 on top of Au FCC(111)
        model = slab(adsorbate=True, surface=surface_facet)
        #view(model)

        # Retrieve index of the C atom
        index = [atom.index for atom in model if atom.symbol == "C"]

        # Mirror model in the x plane with respect to C atom
        # C atom remains in place, the rest of unit cell is shifted accordingly
        model = mirror(model, center_index=index[0], plane='y', surf=surface_facet)
        #view(model)

        # Translate one row of atoms at a time in x and y
        # Move adsorbate to the middle of the unit cell
        model = translation(model, axis=0, surface=surface_facet)
        model = translation(model, axis=1, surface=surface_facet)
        #view(model)


        ### ASSERTION ###
        eps = 1e-8
        # Check if Oxygen moves as expected
        if surface_facet == "100":
            assert (model[19].position - [6.44459559,  0.71187324, 15.07818683] < [eps, eps, eps]).all()
        elif surface_facet == "110":
            assert (model[19].position - [9.48802974,  0.71187324, 14.4695] < [eps, eps, eps]).all()

        ### Rotation example ###
        # TODO: Needs a complete example and assertion - doesn't seem to do anything at first step
        model = rotate_fcc(model, center_index=index[0], surf=surface_facet)
        #view(model)

        ### ASSERTION ###
        # Check if Oxygen moves as expected
        if surface_facet == "100":
            assert (model[19].position - [0.71187324,  8.25040441, 15.07818683] < [eps, eps, eps]).all()
        elif surface_facet == "110":
            assert (model[19].position - [11.29383856,  5.16612676, 14.4695] < [eps, eps, eps]).all()


# Run the example
test_build_symmetry_mirror()
