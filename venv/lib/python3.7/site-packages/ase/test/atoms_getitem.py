from ase.atoms import Atoms
from warnings import warn


w = Atoms('H2O',
          positions=[[2.264, 0.639, 0.876],
                     [0.792, 0.955, 0.608],
                     [1.347, 0.487, 1.234]],
          cell=[3, 3, 3],
          pbc=True)

try:
    print(w[True, False])
except IndexError:
    pass
else:
    # python 3.4 tests skip warnings
    # other python tests are strict
    # warnings will be errors
    warn('')

assert(w[0, 1] == w[True, True, False])
assert(w[0, 1] == w[0:2])

