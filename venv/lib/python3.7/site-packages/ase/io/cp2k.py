# flake8: noqa
"""
Reader for CP2Ks DCD_ALIGNED_CELL format.

Based on [pwtools](https://github.com/elcorto/pwtools).
All information about the dcd format is taken from there.
The way of reading it is also copied from pwtools.
Thanks to Steve for agreeing to this.

Some information also comes directly from the CP2K source,
so if they decide to change anything this here might break.

Some parts are adapted from the extxyz reader.

Contributed by Patrick Melix <chemistry@melix.me>
"""

import numpy as np
from itertools import islice
import os

from ase.atoms import Atoms
from ase.io.formats import index2range

__all__ = ['read_cp2k_dcd', 'iread_cp2k_dcd']

# DCD file header
#   (name, dtype, shape)
# numpy dtypes:
#   i4  = int32
#   f4  = float32 (single precision)
#   f8  = float64 (double precision)
#   S80 = string of length 80 (80 chars)
_HEADER_TYPES = [\
    ('blk0-0',  'i4'      ),  # 84 (start of first block, size=84 bytes)
    ('hdr',     'S4'      ),  # 'CORD'
    ('9int',    ('i4',9)  ),  # 9 ints, mostly 0
    ('timestep','f4'      ),  # timestep (float32)
    ('10int',   ('i4',10) ),  # 10 ints, mostly 0, last is 24
    ('blk0-1',  'i4'      ),  # 84
    ('blk1-0',  'i4'      ),  # 164
    ('ntitle',  'i4'      ),  # 2
    ('remark1', 'S80'     ),  # remark1
    ('remark2', 'S80'     ),  # remark2
    ('blk1-1',  'i4'      ),  # 164
    ('blk2-0',  'i4'      ),  # 4 (4 bytes = int32)
    ('natoms',  'i4'      ),  # natoms (int32)
    ('blk2-1',  'i4'      ),  # 4
    ]

_HEADER_DTYPE = np.dtype(_HEADER_TYPES)

def _bytes_per_timestep(natoms):
    return (4 + 6*8 + 7*4 + 3*4*natoms)

def _read_metainfo(fileobj):
    if not hasattr(fileobj, 'seek'):
        raise TypeError("You should have passed a fileobject opened in binary mode, it seems you did not.")
    fileobj.seek(0)
    header = np.fromfile(fileobj, _HEADER_DTYPE, 1)[0]
    natoms = header['natoms']
    remark1 = str(header['remark1']) #taken from CP2Ks source/motion_utils.F
    if not 'CP2K' in remark1:
        raise ValueError("Header should contain mention of CP2K, are you sure this file was created by CP2K?")

    # dtype for fromfile: nStep times dtype of a timestep data block
    dtype = np.dtype([('x0', 'i4'),
                      ('x1', 'f8', (6,)),
                      ('x2', 'i4', (2,)),
                      ('x3', 'f4', (natoms,)),
                      ('x4', 'i4', (2,)),
                      ('x5', 'f4', (natoms,)),
                      ('x6', 'i4', (2,)),
                      ('x7', 'f4', (natoms,)),
                      ('x8', 'i4')])

    fd_pos = fileobj.tell()
    header_end = fd_pos
    # seek to end
    fileobj.seek(0, os.SEEK_END)
    # number of bytes between fd_pos and end
    fd_rest = fileobj.tell() - fd_pos
    # reset to pos after header
    fileobj.seek(fd_pos)
    # calculate nstep: fd_rest / bytes_per_timestep
    # 4 - initial 48
    # 6*8 - cryst_const_dcd
    # 7*4 - markers between x,y,z and at the end of the block
    # 3*4*natoms - float32 cartesian coords
    nsteps = fd_rest / _bytes_per_timestep(natoms)
    assert fd_rest % _bytes_per_timestep(natoms) == 0 , ("Calculated number of steps is not int, cannot read file")
    nsteps = int(nsteps)
    return dtype, natoms, nsteps, header_end

class DCDChunk:
    def __init__(self, chunk, dtype, natoms, symbols, aligned):
        self.chunk = chunk
        self.dtype = dtype
        self.natoms = natoms
        self.symbols = symbols
        self.aligned = aligned

    def build(self):
        """Convert unprocessed chunk into Atoms."""
        return _read_cp2k_dcd_frame(self.chunk, self.dtype, self.natoms, self.symbols, self.aligned)


def idcdchunks(fd, ref_atoms, aligned):
    """Yield unprocessed chunks for each image."""
    if ref_atoms:
        symbols = ref_atoms.get_chemical_symbols()
    else:
        symbols = None
    dtype, natoms, nsteps, header_end = _read_metainfo(fd)
    bytes_per_step = _bytes_per_timestep(natoms)
    fd.seek(header_end)
    for i in range(nsteps):
        fd.seek(bytes_per_step*i+header_end)
        yield DCDChunk(fd, dtype, natoms, symbols, aligned)


class DCDImageIterator:
    """"""
    def __init__(self, ichunks):
        self.ichunks = ichunks

    def __call__(self, fd, indices=-1, ref_atoms=None, aligned=False):
        self.ref_atoms = ref_atoms
        self.aligned = aligned
        if not hasattr(indices, 'start'):
            if indices < 0:
                indices = slice(indices - 1, indices)
            else:
                indices = slice(indices, indices + 1)

        for chunk in self._getslice(fd, indices):
            yield chunk.build()

    def _getslice(self, fd, indices):
        try:
            iterator = islice(self.ichunks(fd, self.ref_atoms, self.aligned),\
                              indices.start, indices.stop, indices.step)
        except ValueError:
            # Negative indices. Adjust slice to positive numbers.
            dtype, natoms, nsteps, header_end = _read_metainfo(fd)
            indices_tuple = indices.indices(nsteps+1)
            iterator = islice(self.ichunks(fd, self.ref_atoms, self.aligned), *indices_tuple)
        return iterator


iread_cp2k_dcd = DCDImageIterator(idcdchunks)


def read_cp2k_dcd(fileobj, index=-1, ref_atoms=None, aligned=False):
    """ Read a DCD file created by CP2K.

    To yield one Atoms object at a time use ``iread_cp2k_dcd``.

    Note: Other programs use other formats, they are probably not compatible.

    If ref_atoms is not given, all atoms will have symbol 'X'.

    To make sure that this is a dcd file generated with the *DCD_ALIGNED_CELL* key
    in the CP2K input, you need to set ``aligned`` to *True* to get cell information.
    Make sure you do not set it otherwise, the cell will not match the atomic coordinates!
    See (this discussion)[https://groups.google.com/forum/#!searchin/cp2k/Outputting$20cell$20information$20and$20fractional$20coordinates%7Csort:relevance/cp2k/72MhYykrSrQ/5c9Jaw7S9yQJ] for details.
    """
    dtype, natoms, nsteps, header_end = _read_metainfo(fileobj)
    bytes_per_timestep = _bytes_per_timestep(natoms)
    if ref_atoms:
        symbols = ref_atoms.get_chemical_symbols()
    else:
        symbols = ['X' for i in range(natoms)]
    assert natoms == len(symbols), ("Length of ref_atoms does not match natoms from dcd file")
    trbl = index2range(index, nsteps)

    for index in trbl:
        frame_pos = int(header_end + index * bytes_per_timestep)
        fileobj.seek(frame_pos)
        yield _read_cp2k_dcd_frame(fileobj, dtype, natoms, symbols, aligned=aligned)


def _read_cp2k_dcd_frame(fileobj, dtype, natoms, symbols, aligned=False):
    arr = np.fromfile(fileobj, dtype, 1)
    cryst_const = np.empty(6, dtype=np.float64)
    cryst_const[0] = arr['x1'][0,0]
    cryst_const[1] = arr['x1'][0,2]
    cryst_const[2] = arr['x1'][0,5]
    cryst_const[3] = arr['x1'][0,4]
    cryst_const[4] = arr['x1'][0,3]
    cryst_const[5] = arr['x1'][0,1]
    coords = np.empty((natoms,3), dtype=np.float32)
    coords[...,0] = arr['x3'][0,...]
    coords[...,1] = arr['x5'][0,...]
    coords[...,2] = arr['x7'][0,...]
    assert len(coords) == len(symbols)
    if aligned:
        #convention of the cell is (see cp2ks src/particle_methods.F):
        # A is in x
        # B lies in xy plane
        # luckily this is also the ASE convention for Atoms.set_cell()
        atoms = Atoms(symbols, coords, cell=cryst_const, pbc=True)
    else:
        atoms = Atoms(symbols, coords)
    return atoms
