def binary_alloy(model, second_element, n_second_element, random_level=1):
    '''
    Setup a random alloy

    Parameters:

    model: Atoms object
        Atoms object to turn into a binary alloy
    second_element: String
        Name of the new element to introduce
    n_second_element: Integer
        Number of second element species to introduce
    random_level: Integer
        Number of times to randomise element labels (random on random on...)
    '''

    import math, random

    # Prevents unexpected editing of parent object in place
    # Now ensures returned object is different to incoming atoms
    new_model = model.copy()
    count = 0
    labels = new_model.get_chemical_symbols()

    #Using Random Number Generator
    #random.seed(1)

    # Randomly throw in the new element
    while n_second_element > count:
        i = int(math.floor(random.random() * len(new_model)))
        if labels[i] != second_element.capitalize():
            labels[i] = second_element.capitalize()
            count = count + 1

    new_model.set_chemical_symbols(labels)

    # Increase the randomness
    for i in range(random_level):

        count_all = 0
        new_labels = ['']*len(new_model)

        while len(new_model) > count_all:
            j = int(math.floor(random.random()*len(new_model)))
            if new_labels[j] == '':
                new_labels[j] = labels[count_all]
                count_all += 1

        new_model.set_chemical_symbols(new_labels)

    return new_model

def ternary_alloy(model, second_element, third_element, n_second_element, n_third_element, random_level=1):
   
    '''
    Setup a random alloy

    Parameters:

    model: Atoms object
        Atoms object to turn into a binary alloy
    second_element: String
        Name of the new element to introduce
    third_element: String
        Name of the third element in the compounds 
    n_second_element: Integer
        Number of second element species to introduce
    n_third_element: Integer
        Number of third element species to introduce
    random_level: Integer
        Number of times to randomise element labels (random on random on...)
    '''

    import math, random

    # Prevents unexpected editing of parent object in place
    # Now ensures returned object is different to incoming atoms
    new_model = model.copy()
    count = 0
    labels = new_model.get_chemical_symbols()

    #Using Random Number Generator
    #random.seed(1)

   # Randomly throw in the firts new element
    while n_second_element > count:
        i = int(math.floor(random.random() * len(new_model)))
        if labels[i] != second_element.capitalize():
            labels[i] = second_element.capitalize()
            count = count + 1

    new_model.set_chemical_symbols(labels)
    for i in range(random_level):

        count_all = 0
        new_labels = ['']*len(new_model)

        while len(new_model) > count_all:
            j= int(math.floor(random.random()*len(new_model)))
            if new_labels[j] == '':
                new_labels[j] = labels[count_all]
                count_all += 1

       
        
    count = 0
    # Randomly throw the third element to the two elements already presents
    while n_third_element > count:
  
       j = int(math.floor(random.random() * len(new_model)))
       if labels[j] != third_element.capitalize() and labels[j] != second_element.capitalize():
          labels[j] = third_element.capitalize()
          count = count + 1

    new_model.set_chemical_symbols(labels)
    # Increase the randomness
    for i in range(random_level):

        count_all = 0
        new_labels = ['']*len(new_model)

        while len(new_model) > count_all:
            j = int(math.floor(random.random()*len(new_model)))
            if new_labels[j] == '':
                new_labels[j] = labels[count_all]
                count_all += 1

        new_model.set_chemical_symbols(new_labels)

    return new_model


def change_conc(surface, end_conc0, end_conc1, end_conc2, end_conc3, end_conc4, end_conc5, end_conc6):
    '''
      Parameters:

      surface: Atoms object
      end_conc0: Float
           Concentration value for Pd atoms with tag 0
      end_conc1: Float
           Concentration value for Pd atoms with tag 1
      end_conc2: Float
           Concentration value for Pd atoms with tag 2
      end_conc3: Float
           Concentration value for Pd atoms with tag 3
      end_conc4: Float
           Concentration value for Pd atoms with tag 4
      end_conc5: Float
           Concentration value for Pd atoms with tag 5
      end_conc6: Float
           Concentration value for Pd atoms with tag 6
      '''
    import random
    # TODO: make end_conc a list of floats of len(set(atoms.tags))
    end_conc = []
    end_conc.append(end_conc0)
    end_conc.append(end_conc1)
    end_conc.append(end_conc2)
    end_conc.append(end_conc3)
    end_conc.append(end_conc4)
    end_conc.append(end_conc5)
    end_conc.append(end_conc6)

    layer = []
    layer_zn_index = []
    layer_pd_index = []
    layer_no_pd = []
    layer_no_zn = []
    no_loops = []

    # TODO: Specify symbols for substitution
    # TODO: Make the layer count dynamic based on tags
    # TODO: Should the tag assignment be done automatically?
    # TODO: Can work with arrays rather than lists
    # TODO: Include mic and symmetry unique arrangements

    # Calculating the conc. of Pd in layer
    for i in range(0, 7):
        for atom in surface:
            if atom.tag == i:
                layer.append(atom)
                if atom.symbol == "Pd":
                    layer_no_pd.append(atom)
                    layer_pd_index.append(atom.index)
                else:
                    layer_no_zn.append(atom)
                    layer_zn_index.append(atom.index)

        initial_conc = len(layer_no_pd) / len(layer)

        print("The initial concentration of Pd in layer", i, "is:", initial_conc)

        # Changing the conc. to desired conc.
        conc = len(layer_no_pd) / len(layer)

        if conc < end_conc[i]: # Increases conc.
            while conc < end_conc[i]:
                for atom in surface:
                    conc = len(layer_no_pd)/len(layer)
                    if conc < end_conc[i] and atom.symbol == "Zn" and atom.index == random.choice(layer_zn_index):
                        atom.symbol = "Pd"
                        layer_no_pd.append(atom)
                        #print("Atom", atom.index, "is now Pd")
                        conc = len(layer_no_pd) / len(layer)
                        #print("Concentration:", conc)

                no_loops.append("1")
        elif end_conc[i] < conc: # Decreases conc.
            while end_conc[i] < conc:
                for atom in surface:
                    conc = len(layer_no_pd) / len(layer)
                    if end_conc[i] < conc and atom.symbol == "Pd" and atom.index == random.choice(layer_pd_index):
                        atom.symbol = "Zn"
                        del layer_no_pd[0]
                        #print("Atom", atom.index, "is now Pd")
                        conc = len(layer_no_pd) / len(layer)
                        #print("Concentration:", conc)
                no_loops.append("1")

        # Checking calculation
        final_conc = len(layer_no_pd) / len(layer)
        print("The final concentration of Pd in layer", i, "is:", final_conc)
        print("Number of loops:", len(no_loops))
        # Clearing variables
        layer = []
        layer_zn_index = []
        layer_pd_index = []
        layer_no_pd = []
        layer_no_zn = []
        no_loops = []

    # TODO: Come up with sensible return statement
    # TODO: come up with an example for testing