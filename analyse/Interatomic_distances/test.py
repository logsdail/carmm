from ase.data import chemical_symbols, covalent_radii

A = 'H'
B = 'C'

print(covalent_radii[chemical_symbols.index(A)], covalent_radii[chemical_symbols.index(B)])
print(covalent_radii[chemical_symbols.index(A)]+covalent_radii[chemical_symbols.index(B)])