from ase.io.formats import ioformats


text = b"""

  ___ ___ ___ _ _ _  
 |   |   |_  | | | | 
 | | | | | . | | | | 
 |__ |  _|___|_____|  19.8.2b1
 |___|_|             

"""


gpaw = ioformats['gpaw-out']
assert gpaw.match_magic(text)
