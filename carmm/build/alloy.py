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

   # Randomly throw in the first new element
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
   
