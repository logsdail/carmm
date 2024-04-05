"""
Example and test script for the functions seen in hirshfeld.py
These should pull hirshfeld data from an aims.out file
"""


def test_analyse_hirshfeld():

    from carmm.analyse import hirshfeld

    # set the filename, no. of atoms, data name and writing an output file flag
    charge = hirshfeld.extract_hirshfeld(fname='data/mgo/mgo_hirsh_aims.out',
                                         natoms=4,
                                         data='charge',
                                         write=False
                                         )

    assert(charge == (0.34094559, -0.34020775, 0.34094559, -0.34020775))


test_analyse_hirshfeld()



