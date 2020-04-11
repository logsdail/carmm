from ase.data.pubchem import pubchem_search, pubchem_conformer_search
from ase.data.pubchem import pubchem_atoms_search
from ase.data.pubchem import pubchem_atoms_conformer_search

# check class functionality
data = pubchem_search('ammonia', mock_test=True)
atoms = data.get_atoms()
entry_data = data.get_pubchem_data()

# check the various entry styles and the functions that return atoms
atoms = pubchem_search(cid=241, mock_test=True).get_atoms()
atoms = pubchem_atoms_search(smiles='CCOH', mock_test=True)
atoms = pubchem_atoms_conformer_search('octane', mock_test=True)

# check conformer searching
confs = pubchem_conformer_search('octane', mock_test=True)
for conf in confs:
    pass
try:  # check that you can't pass in two args
    atoms = pubchem_search(name='octane', cid=222, mock_test=True)
    raise Exception('Test Failed')
except ValueError:
    pass

try:  # check that you must pass at least one arg
    atoms = pubchem_search(mock_test=True)
    raise Exception('Test Failed')
except ValueError:
    pass
