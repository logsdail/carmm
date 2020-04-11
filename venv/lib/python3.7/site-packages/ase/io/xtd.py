# flake8: noqa
import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom

from ase.utils import basestring
from ase.io.xsd import SetChild, _write_xsd_html
from ase import Atoms

_image_header = ' '*74 + '0.0000\n!DATE     Jan 01 00:00:00 2000\n'
_image_footer = 'end\nend\n'

def _get_atom_str(an,xyz):
    s = '{:<5}'.format(an)
    s += '{:>15.9f}{:>15.9f}{:>15.9f}'.format(xyz[0],xyz[1],xyz[2])
    s += ' XXXX 1      xx      '
    s += '{:<2}'.format(an)
    s += '  0.000\n'
    return s

def write_xtd(filename, images, connectivity=None, moviespeed = 10):
    """Takes Atoms object, and write materials studio file
    atoms: Atoms object
    filename: path of the output file
    moviespeed: speed of animation. between 0 and 10

    note: material studio file cannot use a partial periodic system. If partial
    perodic system was inputted, full periodicity was assumed.
    """
    if moviespeed < 0 or moviespeed > 10:
        raise ValueError('moviespeed only between 0 and 10 allowed')

    if hasattr(images, 'get_positions'):
        images = [images]
        
    XSD,ATR = _write_xsd_html(images,connectivity)
    ATR.attrib['NumChildren'] = '2'
    natoms = len(images[0])
    
    bonds = list()
    if connectivity is not None:
        for i in range(connectivity.shape[0]):
            for j in range(i+1,connectivity.shape[0]):
                if connectivity[i,j]:
                    bonds.append([i,j])

    # non-periodic system
    s = '!BIOSYM archive 3\n'
    if not images[0].pbc.all():
        # Write trajectory
        SetChild(ATR,'Trajectory',{'ID':str(natoms+3+len(bonds)),'Increment':'-1','End':str(len(images)),
            'Type':'arc','Speed':str(moviespeed),'FrameClassType':'Atom'})
        
        # write frame information file
        s += 'PBC=OFF\n'
        for image in images:
            s += _image_header
            s +='\n'
            an = image.get_chemical_symbols()
            xyz = image.get_positions()
            for i in range(natoms):
                s += _get_atom_str(an[i],xyz[i,:])
            s += _image_footer

    # periodic system
    else:
        SetChild(ATR,'Trajectory',{'ID':str(natoms+9+len(bonds)),'Increment':'-1','End':str(len(images)),
            'Type':'arc','Speed':str(moviespeed),'FrameClassType':'Atom'})
            
        # write frame information file
        s += 'PBC=ON\n'
        for image in images:
            s += _image_header
            s += 'PBC'
            vec = image.cell.lengths()
            s+='{:>10.4f}{:>10.4f}{:>10.4f}'.format(vec[0],vec[1],vec[2])
            angles = image.cell.angles() 
            s+='{:>10.4f}{:>10.4f}{:>10.4f}'.format(angles[0],angles[1],angles[2])
            s+='\n'
            an = image.get_chemical_symbols()

            angrad = np.deg2rad(angles)
            cell = np.zeros((3,3))
            cell[0,:] = [vec[0],0 ,0]
            cell[1,:] = np.array([np.cos(angrad[2]), np.sin(angrad[2]), 0]) * vec[1]
            cell[2,0] = vec[2]*np.cos(angrad[1])
            cell[2,1] = (vec[1]*vec[2]*np.cos(angrad[0])-cell[1,0]*cell[2,0])/cell[1,1]
            cell[2,2] = np.sqrt(vec[2]**2-cell[2,0]**2-cell[2,1]**2)
            xyz = np.dot(image.get_scaled_positions(),cell)
            for i in range(natoms):
                s += _get_atom_str(an[i],xyz[i,:])
            s += _image_footer

    # print arc file
    if isinstance(filename,str):
        farcname = filename[:-3] + 'arc'
    else:
        farcname = filename.name[:-3] + 'arc'
    farc = open(farcname, 'w')
    farc.write(s)
    farc.close()

    # check if file is an object or not.
    if isinstance(filename, basestring):
        f = open(filename, 'w')
    else:  # Assume it's a 'file-like object'
        f = filename

    # Return a pretty-printed XML string for the Element.
    rough_string = ET.tostring(XSD, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    Document = reparsed.toprettyxml(indent='\t')

    # write
    f.write(Document)
    f.close()
    
def read_xtd(filename, index=-1):
    """Import xtd file (Materials Studio)
    
    Xtd files always come with arc file, and arc file
    contains all the relevant information to make atoms
    so only Arc file needs to be read
    """
    if isinstance(filename,str):
        f = filename[:-3] + 'arc'
    else:
        f = filename.name[:-3] + 'arc'
    f = open(f,'r')
    images = list()

    # the first line is comment
    f.readline()
    pbc = 'ON' in f.readline()
    l = f.readline()
    while l != '':
        if '!' not in l: # flag for the start of an image
            l = f.readline()
            continue
        if pbc:
            l = f.readline()
            cell = [float(d) for d in l.split()[1:]]
        else:
            f.readline()
        symbols = []
        coords = []
        while True:
            l = f.readline().split()
            if 'end' in l:
                break
            symbols.append(l[0])
            coords.append([float(x) for x in l[1:4]])
        if pbc:
            image = Atoms(symbols, positions=coords, cell=cell, pbc=pbc)
        else:
            image = Atoms(symbols, positions=coords, pbc=pbc)
        images.append(image)
        l = f.readline()
        
                
    if not index:
        return images
    else:
        return images[index]

    

