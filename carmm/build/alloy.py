def binary_alloy(atoms, second_element, n_second_element, random_level=1):
    '''
    Setup a random alloy

    Parameters:

    atoms: Atoms object
        Atoms object to turn into a binary alloy
    second_element: String
        Name of the new element to introduce
    n_second_element: Integer
        Number of second element species to introduce
    random_level: Integer
        Number of times to randomise element labels (random on random on...)
    '''

    import math, random

    count = 0
    labels = atoms.get_chemical_symbols()

    #Using Random Number Generator
    #random.seed(1)

    # Randomly throw in the new element
    while n_second_element > count:
        i = int(math.floor(random.random() * len(atoms)))
        if labels[i] != second_element.capitalize():
            labels[i] = second_element.capitalize()
            count = count + 1

    atoms.set_chemical_symbols(labels)

    # Increase the randomness
    for i in range(random_level):

        count_all = 0
        new_labels = ['']*len(atoms)

        while len(atoms) > count_all:
            j = int(math.floor(random.random()*len(atoms)))
            if new_labels[j] == '':
                new_labels[j] = labels[count_all]
                count_all += 1

        atoms.set_chemical_symbols(new_labels)

    return atoms

def ternary_alloy(atoms, second_element, third_element, n_second_element, n_third_element, random_level=1):
   
    '''
    Setup a random alloy

    Parameters:

    atoms: Atoms object
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

    count = 0
    labels = atoms.get_chemical_symbols()

    #Using Random Number Generator
    #random.seed(1)

   # Randomly throw in the firts new element
    while n_second_element > count:
        i = int(math.floor(random.random() * len(atoms)))
        if labels[i] != second_element.capitalize():
            labels[i] = second_element.capitalize()
            count = count + 1

    atoms.set_chemical_symbols(labels)
    for i in range(random_level):

        count_all = 0
        new_labels = ['']*len(atoms)

        while len(atoms) > count_all:
            j= int(math.floor(random.random()*len(atoms)))
            if new_labels[j] == '':
                new_labels[j] = labels[count_all]
                count_all += 1

       
        
    count = 0
    # Randomly throw the third element to the two elements already presents
    while n_third_element > count:
  
       j = int(math.floor(random.random() * len(atoms)))
       if labels[j] != third_element.capitalize():
          labels[j] = third_element.capitalize()
          count = count + 1

    atoms.set_chemical_symbols(labels)
    # Increase the randomness
    for i in range(random_level):

        count_all = 0
        new_labels = ['']*len(atoms)

        while len(atoms) > count_all:
            j = int(math.floor(random.random()*len(atoms)))
            if new_labels[j] == '':
                new_labels[j] = labels[count_all]
                count_all += 1

        atoms.set_chemical_symbols(new_labels)

    return atoms
   
