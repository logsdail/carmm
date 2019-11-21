import sys
from gpaw import GPAW
from gpaw import restart
from ase.io import write
from gpaw.elf import ELF

filename = sys.argv[1]
atoms, calc = restart(filename)
calc.set_positions()

elf = ELF(calc)
#elf.initialize(calc)
#elf.update(calc.wfs)
elf.update()

elf_out = elf.get_electronic_localization_function(gridrefinement=2)

## Write to file
## write(filename+'.plt', atoms, data=elf_out)
write(filename+'.cube', atoms, data=elf_out)
