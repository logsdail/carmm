from ase.lattice import bravais_classes


_lattice_system = ('Øaammmmmmmmmmmmmoooooooooooooooooooooooooooooooooooooooooo'
                   'ooooooooooooooooottttttttttttttttttttttttttttttttttttttttt'
                   'ttttttttttttttttttttttttttthhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh'
                   'hhhhhhhhhhhhhhhhhhhhhcccccccccccccccccccccccccccccccccccc')

_lattice_centering = ('ØPPPPCPPCCPPCPPCPPPPCCFIIPPPPPPPPPPCCCCCCCFFIIIPPPPPPPP'
                      'PPPPPPPPCCCCCCFFIIIIPPPPIIPIPPPPIIPPPPPPPPIIPPPPPPPPII'
                      'IIPPPPPPPPIIIIPPPPPPPPPPPPPPPPIIIIPPPRPRPPPPPPRPPPPRRP'
                      'PPPRRPPPPPPPPPPPPPPPPPPPPPPPPPPPPFIPIPPFFIPIPPFFIPPIPF'
                      'IPFIPPPPFFFFII')


def get_bravais_class(sg):

    sg = int(sg)
    if sg < 1:
        raise ValueError('Spacegroup must be positive, but is {}'.format(sg))
    if sg > 230:
        raise ValueError('Bad spacegroup', sg)
    pearson_symbol = _lattice_system[sg] + _lattice_centering[sg]
    return bravais_classes[pearson_symbol]
