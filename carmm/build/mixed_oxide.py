from ase import Atoms

def mixed_oxide(atoms_object:Atoms, metal1:str, metal2:str, ratio:int, metal1_charge:int, metal2_charge:int, 
                metal1_moment:int, metal2_moment:int, up_or_down_spin=None):
    '''
    Provided an Atoms object of a metal oxide with formula AxOy, this function creates another
    Atoms object for a mixed oxide with formula Ax-kBkOy where A and B are the two metals.
    x is the initial stoichiometry of A in AxOy.
    The value of k can be figured out from the equation ratio(metal_A:metal_B) = x-k/k for a given value of ratio
    IMPORTANT NOTE: As of now this functionality only works for oxides which the metal A is in only one oxidation state
    like CoO or Mn2O3, etc.

    :param atoms_object: the Atoms object for the oxide AxOy with initial moments and charges specified. This could be obtained by reading a geometry
    file.
    :param metal1: (String) The metal A. eg: Co, Mn, Fe, etc
    :param metal2: (String) The metal B in Ax-kBkOy
    :param ratio: (int) The desired ratio of metal A to metal B in the output geometry. In this case x-k/k
    :param metal1_charge: (int) The charge on metal A that you have specified in the input geometry.in
    :param metal2_charge: (int) The charge on metal B that you want in the output geometry_mixed.in.
    NOTE: make sure that both metal A and B both have same charges (as this may result in a charge imbalance in the
    final structure
    :param metal1_moment: (int) The initial moment on metal A that you have specified in the input geometry.in. No need to specify the direction
    :param metal2_moment: (int) The initial moment on metal B that you desire to have in the mixed oxide geometry. No need to specify the direction
    :param up_or_down_spin: if the metal1 has zero magnetic moment and metal2 has a magnetic moment specified by
    metal2_moment, please specify the sign for the magnetic moment of the metal2. Positive (up-spin) or Negative (down-spin)
    :return: Returns an atoms object of the mixed oxide geometry and a geometry.in file could be written using this atoms object
    '''
    
    import math
    import numpy as np
    from ase import Atoms
    import random

    if up_or_down_spin is not None:
        up_or_down_spin=up_or_down_spin.capitalize()
        print(up_or_down_spin)

    atoms = atoms_object
    sym = np.array(atoms.get_chemical_symbols())
    charge= np.array(atoms.get_initial_charges())
    mom = np.array(atoms.get_initial_magnetic_moments())
    mom_sign= np.sign(mom)
    m1_ind = []
    for ind,symbol in enumerate(sym):
        if symbol==metal1.capitalize() and charge[ind]==metal1_charge and np.abs(mom[ind])==metal1_moment:
            m1_ind.append(ind)

    # determining the number of metal1 ions to be replaced with metal2 which depends on the ratio provided by the user.
    k = math.floor((len(m1_ind)/(1+ratio)))
    # random replacement of metal1 ion with metal2 ions
    m2_ind = random.sample(m1_ind, k)
    '''
    This functionality relies on the random replacement of metal 1 ions with metal 2 ions. Strictly speaking, randomly 
    generated geometry files for a given ratio must be tested for energy and it should be figured out which one is the 
    most stable. 
    But this combinatorial problem can be very huge and hence expensive.'''
    for i in m2_ind:
        sym[i] = metal2.capitalize()
        charge[i] = metal2_charge
        if mom_sign[ind]==0:
            if up_or_down_spin==None:
                raise Exception('the magnetic moment of the metal atom that you desire to replace is zero'
                                '(non-magnetic). Please provide the value for the argument up_or_down_spin')
            elif up_or_down_spin=='Positive':
                print('the magnetic moment of the metal atom that you desire to replace is zero(non-magnetic). '
                      'Using a positive (up-spin) moment')
                mom[i] = metal2_moment
            elif up_or_down_spin == 'Negative':
                print('the magnetic moment of the metal atom that you desire to replace is zero(non-magnetic). '
                      'Using a negative (down-spin) moment')
                mom[i] = (-1)*metal2_moment
            else:
                raise Exception(f'Invalid input to the argument up_or_down_spin. Should be either \'positive\' or '
                                f'\'negative\'. You provided up_or_down_spin={up_or_down_spin}')
        elif mom_sign[ind]!=0:
            mom[i] = mom_sign[ind]*metal2_moment

    atoms.set_chemical_symbols(sym)
    atoms.set_initial_charges(charge)
    atoms.set_initial_magnetic_moments(mom)

    if neutrality_check(atoms):
        mixed_oxide = atoms
        print('The mixed oxide geometry has been created as an atoms object and a geometry.in file can be written'
              'using the write functionality of ASE.')
        return mixed_oxide
    elif not neutrality_check(atoms):
        raise Exception(f'The mixed oxide generated is not neutral and has an excess charge.'
                        f'Please make sure that the metal1 and metal1 have the same charge')

def neutrality_check(atoms_object: Atoms):
    import numpy as np
    charge_list = np.array(atoms_object.get_initial_charges())
    charge_sum = np.sum(charge_list)
    if charge_sum == 0:
        print('Neutrality check passed')
        return True
    else:
        return False
