'''
The code allows to calculate the energy and forces of two fragments within the same system
with two different calculators. The idea is to use mace calculator for pre-optimisation in surface 
reactions, where adsorbates can be represented by a seperate foundation model than the rest of the system. 
The same script may be used for supported nanocluster adsorption/reaction in a 3 calculator setup.

TODO:
    1. Integrate with carmm workflow
    2. Generalise the code to take system and define fragments
    3. Include the fragments to contain some nearest neighbour indices from other fragments for
       smooth distribution of forces/energy
    4. validate stress parameter
'''

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator

class PartitionedCalculator(Calculator):
    def __init__(self, calculators, atom_groups, **kwargs):
        """Calculator that partitions atoms into groups with different calculators
        
        Args:
            calculators (list): List of ASE calculator instances
            atom_groups (list): List of atom index lists for each calculator
        """
        super().__init__(**kwargs)
        self.calculators = calculators
        self.atom_groups = atom_groups
        
    def calculate(self, atoms=None, properties=None, system_changes=None):
        super().calculate(atoms, properties, system_changes)
        
        properties = properties or []
        total_energy = 0.0
        forces = np.zeros((len(atoms), 3))
        stresses = np.zeros(6) if "stress" in properties else None
        
        for calc, indices in zip(self.calculators, self.atom_groups):
            # Create subsystem for this calculator
            subsystem = atoms[indices]
            subsystem.calc = calc
            
            # Calculate properties
            energy = subsystem.get_potential_energy()
            force = subsystem.get_forces()
            
            # Aggregate results
            total_energy += energy
            forces[indices] += force
            if 'stress' in properties:
                stress = subsystem.get_stress()
                stresses += stress * (len(indices) / len(atoms))  # Approximate stress partitioning
                
        self.results = {
            'energy': total_energy,
            'forces': forces,
        }
        if "stress" in properties:
            self.results["stress"] = stresses        
        return self.results

# Example usage:
from ase.calculators.emt import EMT
#from ase.calculators.lj import LennardJones
#from mace.calculators import mace_mp

#mp_small = mace_mp(model="small", default_dtype="float64", device='cpu')

# Create a molecule with 4 atoms
atoms = Atoms('H2O2', positions=[[0, 0, 0], [1, 0, 0], [0.5, 1, 0], [0.5, 0, 1]])
atoms1 =atoms.copy()

#Define the fragments within the molecule
frag1=[0,1]
frag2=[2,3]

# Define calculator groups (first 2 atoms with EMT, rest with EMT again)
calc = PartitionedCalculator(
    calculators=[EMT(), EMT()],      #swap one of the calc with for different calculators
    atom_groups=[frag1, frag2]
)

# Perform single-point calculation
props = calc.calculate(atoms)
print(props['forces'])

#check the result against a single calculator (EMT in this case)
atoms1.calc=EMT()
print(atoms1.get_forces())