from ase.visualize import view
from data.model_gen import get_example_slab as slab
from carmm.analyse.neb_tools.pathway import minimize_distance, apply_sequence
import numpy as np


def test_neb_path():
    '''
    A function showcasing the use of NEB pathway generator on FCC low-index surfaces. The shift in atomic coordinates
    is minimised based on symmetry and periodicity, resulting in a simpler pathway, less complex potential energy surface,
    and overall easier calculation. Supported facets: "111", "100", "110". Assumption - single sided adsorption.

    '''

    for surface_facet in ["111", "100", "110"]:
        # Toy model of CO2 on top of Au FCC(111)
        initial = slab(adsorbate=True, surface=surface_facet)
        final = slab(adsorbate=True, surface=surface_facet)

        for n in [1,2,3]:
            initial[-n].position += (-4, -4, 0)

        list_of_geometries = minimize_distance(initial, final, (3,3,2), surface_facet)

        # initial, final = apply_sequence(initial, final, list_of_geometries[0][0], surface_facet)
        # view([initial, final])


        assertion_values = {"111_seq":[2,1,0], "111_d":15.280125,
                            "100_seq":[1,2,0], "100_d":12.311127,
                            "110_seq":[1,0,0], "110_d":14.564399}

        assert (np.array(assertion_values[surface_facet+"_seq"]) == list_of_geometries[0][0]).all()
        assert assertion_values[surface_facet+"_d"] - list_of_geometries[0][1] < 1e-6


test_neb_path()

