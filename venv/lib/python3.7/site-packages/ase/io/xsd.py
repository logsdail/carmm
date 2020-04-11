# flake8: noqa
import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom

from ase.utils import basestring
from ase import Atoms


def read_xsd(fd):
    tree = ET.parse(fd)
    root = tree.getroot()

    atomtreeroot = root.find('AtomisticTreeRoot')
    # if periodic system
    if atomtreeroot.find('SymmetrySystem') is not None:
        symmetrysystem = atomtreeroot.find('SymmetrySystem')
        mappingset = symmetrysystem.find('MappingSet')
        mappingfamily = mappingset.find('MappingFamily')
        system = mappingfamily.find('IdentityMapping')

        coords = list()
        cell = list()
        formula = str()

        for atom in system:
            if atom.tag == 'Atom3d':
                symbol = atom.get('Components')
                formula += symbol

                xyz = atom.get('XYZ')
                if xyz:
                    coord = [float(coord) for coord in xyz.split(',')]
                else:
                    coord = [0.0, 0.0, 0.0]
                coords.append(coord)
            elif atom.tag == 'SpaceGroup':
                avec = [float(vec) for vec in atom.get('AVector').split(',')]
                bvec = [float(vec) for vec in atom.get('BVector').split(',')]
                cvec = [float(vec) for vec in atom.get('CVector').split(',')]

                cell.append(avec)
                cell.append(bvec)
                cell.append(cvec)

        atoms = Atoms(formula, cell=cell, pbc=True)
        atoms.set_scaled_positions(coords)
        return atoms
        # if non-periodic system
    elif atomtreeroot.find('Molecule') is not None:
        system = atomtreeroot.find('Molecule')

        coords = list()
        formula = str()

        for atom in system:
            if atom.tag == 'Atom3d':
                symbol = atom.get('Components')
                formula += symbol

                xyz = atom.get('XYZ')
                coord = [float(coord) for coord in xyz.split(',')]
                coords.append(coord)

        atoms = Atoms(formula, pbc=False)
        atoms.set_scaled_positions(coords)
        return atoms


def CPK_or_BnS(element):
    """Determine how atom is visualized"""
    if element in ['C', 'H', 'O', 'S', 'N']:
        visualization_choice = 'Ball and Stick'
    else:
        visualization_choice = 'CPK'
    return visualization_choice

def SetChild(parent,childname,props):
    Child = ET.SubElement(parent,childname)
    for key in props:
        Child.set(key,props[key])
    return Child

def SetBasicChilds():
    """
    Basic property setup for Material Studio File
    """
    XSD = ET.Element('XSD')
    XSD.set('Version', '6.0')

    ATR = SetChild(XSD,'AtomisticTreeRoot',{'ID':'1','NumProperties':'40','NumChildren':'1'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'AngleEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'BendBendEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'BendTorsionBendEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'BondEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'Atom','Name':'EFGAsymmetry','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'Atom','Name':'EFGQuadrupolarCoupling','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'ElectrostaticEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'GrowthFace','Name':'FaceMillerIndex','Type':'MillerIndex'})
    SetChild(ATR,'Property',{'DefinedOn':'GrowthFace','Name':'FacetTransparency','Type':'Float'})
    SetChild(ATR,'Property',{'DefinedOn':'Bondable','Name':'Force','Type':'CoDirection'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'HydrogenBondEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'Bondable','Name':'ImportOrder','Type':'UnsignedInteger'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'InversionEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'Atom','Name':'IsBackboneAtom','Type':'Boolean'})
    SetChild(ATR,'Property',{'DefinedOn':'Atom','Name':'IsChiralCenter','Type':'Boolean'})
    SetChild(ATR,'Property',{'DefinedOn':'Atom','Name':'IsOutOfPlane','Type':'Boolean'})
    SetChild(ATR,'Property',{'DefinedOn':'BestFitLineMonitor','Name':'LineExtentPadding','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'Linkage','Name':'LinkageGroupName','Type':'String'})
    SetChild(ATR,'Property',{'DefinedOn':'PropertyList','Name':'ListIdentifier','Type':'String'})
    SetChild(ATR,'Property',{'DefinedOn':'Atom','Name':'NMRShielding','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'NonBondEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'Bondable','Name':'NormalMode','Type':'Direction'})
    SetChild(ATR,'Property',{'DefinedOn':'Bondable','Name':'NormalModeFrequency','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'Bondable','Name':'OrbitalCutoffRadius','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'BestFitPlaneMonitor','Name':'PlaneExtentPadding','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'PotentialEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ScalarFieldBase','Name':'QuantizationValue','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'RestraintEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'SeparatedStretchStretchEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'Trajectory','Name':'SimulationStep','Type':'Integer'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'StretchBendStretchEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'StretchStretchEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'StretchTorsionStretchEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'TorsionBendBendEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'TorsionEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'TorsionStretchEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'ValenceCrossTermEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'ValenceDiagonalEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'ClassicalEnergyHolder','Name':'VanDerWaalsEnergy','Type':'Double'})
    SetChild(ATR,'Property',{'DefinedOn':'SymmetrySystem','Name':'_Stress','Type':'Matrix'})
    return ATR, XSD

def _write_xsd_html(images,connectivity=None):
    ATR, XSD = SetBasicChilds()
    natoms = len(images[0])
    atom_element = images[0].get_chemical_symbols()
    atom_cell = images[0].get_cell()
    atom_positions = images[0].get_positions()
    # Set up bonds
    bonds = list()
    if connectivity is not None:
        for i in range(connectivity.shape[0]):
            for j in range(i+1,connectivity.shape[0]):
                if connectivity[i,j]:
                    bonds.append([i,j])

    # non-periodic system
    if not images[0].pbc.all():
        Molecule = SetChild(ATR,'Molecule',{'ID':'2','NumChildren':str(natoms+len(bonds)),'Name':'Lattice=&quot1.0'})
        # writing images[0]
        for x in range(natoms):
            Props = {'ID': str(x + 3),'Name':(atom_element[x] + str(x + 1)),'UserID':str(x + 1),'DisplayStyle':CPK_or_BnS(atom_element[x])}
            Props['XYZ']='%1.16f,%1.16f,%1.16f' %(atom_positions[x,0],atom_positions[x,1],atom_positions[x,2])
            Props['Components']=atom_element[x]
            bondstr = ['%i'% (i + 3 + natoms) for i,bond in enumerate(bonds) if x in bond]
            if bondstr:
                Props['Connections']=','.join(bondstr)
            SetChild(Molecule,'Atom3d',Props)
        for x in range(len(bonds)):
            SetChild(Molecule,'Bond',{'ID':str(x + 3 + natoms),'Connects':'%i,%i'%(bonds[x][0] + 3,bonds[x][1] + 3)})
    # periodic system
    else:
        atom_positions = np.dot(atom_positions, np.linalg.inv(atom_cell))
        Props = {}
        Props['ID']='2'
        Props['Mapping']= '3'
        Props['Children']= ''.join(['%1.0f,' % (x) for x in range(4, natoms + len(bonds) + 4)]+[str(natoms + len(bonds) + 4)])
        Props['Normalized']= '1'
        Props['Name']= 'SymmSys'
        Props['UserID']= str(natoms + 18)
        Props['XYZ']= '0.00000000000000,0.00000000000000,0.000000000000000'
        Props['OverspecificationTolerance']= '0.05'
        Props['PeriodicDisplayType']= 'Original'
        SymmSys = SetChild(ATR,'SymmetrySystem',Props)

        Props = {}
        Props['ID']= str(natoms + len(bonds) + 5)
        Props['SymmetryDefinition']= str(natoms + 4)
        Props['ActiveSystem']= '2'
        Props['NumFamilies']= '1'
        Props['OwnsTotalConstraintMapping']= '1'
        Props['TotalConstraintMapping']= '3'
        MappngSet = SetChild(SymmSys, 'MappingSet',Props)

        Props = {}
        Props['ID']= str(natoms + len(bonds) + 6)
        Props['NumImageMappings']= '0'
        MappngFamily = SetChild(MappngSet, 'MappingFamily',Props)

        Props = {}
        Props['ID']= str(natoms + len(bonds) + 7)
        Props['Element']= '1,0,0,0,0,1,0,0,0,0,1,0'
        Props['Constraint']= '1,0,0,0,0,1,0,0,0,0,1,0'
        Props['MappedObjects']= ','.join(['%1.0f' % (x) for x in range(4, natoms + len(bonds) + 4)])
        Props['DefectObjects']= str(natoms + len(bonds) + 4) + ',' + str(natoms + len(bonds) + 8)
        Props['NumImages']= str(natoms + len(bonds))
        Props['NumDefects']= '2'
        IdentMappng = SetChild(MappngFamily, 'IdentityMapping',Props)

        SetChild(MappngFamily,'MappingRepairs',{'NumRepairs':'0'})

        # writing atoms
        for x in range(natoms):
            Props = {}
            Props['ID']= str(x + 4)
            Props['Mapping']=  str(natoms + len(bonds) + 7)
            Props['Parent']= '2'
            Props['Name']= (atom_element[x] + str(x + 1))
            Props['UserID']= str(x + 1)
            Props['DisplayStyle']= CPK_or_BnS(atom_element[x])
            Props['Components']= atom_element[x]
            Props['XYZ']='%1.16f,%1.16f,%1.16f' %(atom_positions[x,0],atom_positions[x,1],atom_positions[x,2])
            bondstr = ['%i'% (i + 4 + natoms + 1) for i,bond in enumerate(bonds) if x in bond]
            if bondstr:
                Props['Connections']=','.join(bondstr)
            SetChild(IdentMappng, 'Atom3d',Props)

        for x in range(len(bonds)):
            SetChild(IdentMappng,'Bond',{'ID':str(x + 4 + natoms + 1),'Mapping':str(natoms + len(bonds) + 7),
                'Parent':'2','Connects':'%i,%i'%(bonds[x][0] + 4,bonds[x][1] + 4)})

        Props={}
        Props['ID']= str(natoms + 4)
        Props['Parent']= '2'
        Props['Children']= str(natoms + len(bonds) + 8)
        Props['DisplayStyle']= 'Solid'
        Props['XYZ']= '0.00,0.00,0.00'
        Props['Color']= '0,0,0,0'
        Props['AVector']= ','.join(['%1.16f' % atom_cell[0, x] for x in range(3)])
        Props['BVector']= ','.join(['%1.16f' % atom_cell[1, x] for x in range(3)])
        Props['CVector']= ','.join(['%1.16f' % atom_cell[2, x] for x in range(3)])
        Props['OrientationBase']= 'C along Z, B in YZ plane'
        Props['Centering']= '3D Primitive-Centered'
        Props['Lattice']= '3D Triclinic'
        Props['GroupName']= 'GroupName'
        Props['Operators']= '1,0,0,0,0,1,0,0,0,0,1,0'
        Props['DisplayRange']= '0,1,0,1,0,1'
        Props['LineThickness']= '2'
        Props['CylinderRadius']= '0.2'
        Props['LabelAxes']= '1'
        Props['ActiveSystem']= '2'
        Props['ITNumber']= '1'
        Props['LongName']= 'P 1'
        Props['Qualifier']= 'Origin-1'
        Props['SchoenfliesName']= 'C1-1'
        Props['System']= 'Triclinic'
        Props['Class']= '1'
        SetChild(IdentMappng, 'SpaceGroup',Props)

        SetChild(IdentMappng, 'ReciprocalLattice3D',{'ID':str(natoms + len(bonds) + 8),
            'Parent':str(natoms + 4)})

        SetChild(MappngSet, 'InfiniteMapping',{'ID':'3','Element':'1,0,0,0,0,1,0,0,0,0,1,0',
            'MappedObjects':'2'})

    return XSD,ATR

def write_xsd(filename, images, connectivity=None):
    """Takes Atoms object, and write materials studio file
    atoms: Atoms object
    filename: path of the output file
    connectivity: number of atoms by number of atoms matrix for connectivity
    between atoms (0 not connected, 1 connected)

    note: material studio file cannot use a partial periodic system. If partial
    perodic system was inputted, full periodicity was assumed.
    """
    if hasattr(images, 'get_positions'):
        images = [images]

    XSD,ATR = _write_xsd_html(images,connectivity)

    # check if file is an object or not.
    if isinstance(filename, basestring):
        f = open(filename, 'w')
    else:  # Assume it's a 'file-like object'
        f = filename

    # Return a pretty-printed XML string for the Element.
    rough_string = ET.tostring(XSD, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    Document = reparsed.toprettyxml(indent='\t')

    f.write(Document)
