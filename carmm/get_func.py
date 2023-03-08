import os
import sys



def list_functions(folder_path, output_file_path, relative_path=''):
    """
    Recursively list all functions in modules relative to the given folder path
    and write the function names into the given output file path.
    """
    # Add the top-level folder path to the module search path
    top_level_folder = os.path.abspath(os.path.join(folder_path, '..', '..'))
    sys.path.insert(0, top_level_folder)

    # Open the output file in append mode
    with open(output_file_path, 'a') as output_file:
        # Iterate over all files and subdirectories in the given folder path
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            # If the file is a directory, recursively call this function on it
            if os.path.isdir(file_path):
                new_relative_path = os.path.join(relative_path, filename)
                list_functions(file_path, output_file_path, new_relative_path)
            # If the file is a Python module, inspect it for functions
            elif filename.endswith('.py'):
                # Get the module name without the .py extension
                module_name = os.path.splitext(filename)[0]
                # Import the module
                module_path = os.path.join(relative_path, module_name).replace(os.sep, '.')
                module = __import__(module_path, fromlist=[''])
                # Get the names of all objects defined in the module
                object_names = dir(module)
                # Write the names of functions to the output file
                for object_name in object_names:
                    object = getattr(module, object_name)
                    if callable(object) and hasattr(object, '__module__') and object.__module__ == module_path:
                        # output_file.write(f'{os.path.join(relative_path, module_name)}.{object_name}\n')
                        # Add a double space at the end of the line for .md markdown
                        function = ".".join(relative_path.split("\\") + [module_name, object_name + '  \n'])
                        set_of_functions.add(function)


    # Remove the top-level folder path from the module search path
    sys.path.pop(0)


set_of_functions = set()
filename ='FUNCTIONS.md'
list_functions('.', filename)
list_of_functions = sorted(list(set_of_functions))

with open(filename, 'a') as file:

    [file.write(line + "  ") for line in list_of_functions]
