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
                                         write=True,
                                         outname='data/mgo/hirshfeld.txt'
                                         )

    assert(charge == [0.34094559, -0.34020775, 0.34094559, -0.34020775])

    vol = hirshfeld.extract_hirshfeld(fname='data/mgo/mgo_hirsh_aims.out',
                                         natoms=4,
                                         data='volume',
                                         write=False
                                         )

    assert(vol == [84.98662789, 22.60472263, 84.98662789, 22.60472263])

    volf = hirshfeld.extract_hirshfeld(fname='data/mgo/mgo_hirsh_aims.out',
                                         natoms=4,
                                         data='volume f',
                                         write=False
                                         )

    assert(volf == [93.92897117, 23.60491416, 93.92897117, 23.60491416])

    dipv = hirshfeld.extract_hirshfeld(fname='data/mgo/mgo_hirsh_aims.out',
                                         natoms=4,
                                         data='dipole vector',
                                         write=False
                                         )

    assert(dipv == [[0.0, 0.0, -0.0], [0.0, 0.0, -0.0], [-0.0, -0.0, 0.0], [-0.0, -0.0, 0.0]])

    dipm = hirshfeld.extract_hirshfeld(fname='data/mgo/mgo_hirsh_aims.out',
                                         natoms=4,
                                         data='dipole moment',
                                         write=False
                                         )

    assert(dipm == [0.0, 0.0, 0.0, 0.0])

    secn = hirshfeld.extract_hirshfeld(fname='data/mgo/mgo_hirsh_aims.out',
                                         natoms=4,
                                         data='second',
                                         write=False
                                         )

    assert(secn == [[0.22597625, 0.0, 0.0, 0.0, 0.22597621, 0.0, 0.0, 0.0, 0.2249063],
                    [-0.00625081, 0.0, -0.0, 0.0, -0.00625081, -0.0, -0.0, -0.0, -0.00618102],
                    [0.22597625, 0.0, 0.0, 0.0, 0.22597621, 0.0, 0.0, 0.0, 0.2249063],
                    [-0.00625081, 0.0, -0.0, 0.0, -0.00625081, -0.0, -0.0, -0.0, -0.00618102]])

    # Test for VMD output

    import numpy as np

    hirshfeld.vmd_out(np.array(charge), fname='data/mgo/vmd_chrgs.txt')
    with open('data/mgo/vmd_chrgs.txt','r') as vmd:
        lines = vmd.readlines()
    assert(lines == ['0.34094559\n', '-0.34020775\n', '0.34094559\n', '-0.34020775\n'])

test_analyse_hirshfeld()



