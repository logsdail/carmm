"""Reads chemical data in MDL Molfile format. 

See https://en.wikipedia.org/wiki/Chemical_table_file
"""
from ase.atoms import Atoms


def read_mol(fileobj):
    lines = fileobj.readlines()
    L1 = lines[3]

    # The V2000 dialect uses a fixed field length of 3, which means there
    # won't be space between the numbers if there are 100+ atoms, and
    # the format doesn't support 1000+ atoms at all.
    if L1.rstrip().endswith('V2000'):
        natoms = int(L1[:3].strip())
    else:
        natoms = int(L1.split()[0])
    positions = []
    symbols = []
    for line in lines[4:4 + natoms]:
        x, y, z, symbol = line.split()[:4]
        symbols.append(symbol)
        positions.append([float(x), float(y), float(z)])
    return Atoms(symbols=symbols, positions=positions)
