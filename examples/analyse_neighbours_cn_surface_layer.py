def test_analyse_neighbours_cn_surface_layer():

    from carmm.analyse.neighbours import cn_surface_layers
    from examples.data.model_gen import get_example_slab

    # Get example fcc111 Au slab with 2 Cu atoms on top
    model = get_example_slab(adsorbate=True, type="2Cu")

    # Ensure Cu atoms are close to the surface
    model[-1].z -= 1
    model[-2].z -= 1

    # set verbose to True to see printed data
    cn_per_atom, cn_per_layer = cn_surface_layers(model, verbose=False)

    # Assert that dynamically assigned symbols are working
    assert [i for i in cn_per_atom[0]] == ['symbol', 'index', 'layer', 'Au_neighbors', 'Cu_neighbors']
    assert [i for i in cn_per_layer[0]] == ['layer', 'Au_neighboring_w_Au', 'Au_neighboring_w_Cu', 'Cu_neighboring_w_Au',
                                            'Cu_neighboring_w_Cu', 'Au_concentration_per_layer', 'Cu_concentration_per_layer']
    assert len(cn_per_layer[0]) == 7

    # Assert that extracted values are correct
    assert cn_per_layer[0]["Au_neighboring_w_Cu"] == 3.0
    assert cn_per_layer[1]["Au_neighboring_w_Au"] == 9.0

test_analyse_neighbours_cn_surface_layer()
