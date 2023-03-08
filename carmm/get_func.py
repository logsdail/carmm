import os
import importlib

import os
import importlib

def list_functions_in_folder(folder_path, output_file):
    # Loop over all the files in the folder
    for filename in os.listdir(folder_path):
        # Get the full path of the file
        file_path = os.path.join(folder_path, filename)
        # Check if the file is a Python module
        if filename.endswith('.py'):
            # Remove the '.py' extension to get the module name
            module_name = filename[:-3]
            # Prepend the module name with the full path
            module_path = os.path.relpath(folder_path)
            module_name_with_path = os.path.join(module_path, module_name).replace('\\', '/')
            # Import the module
            module = importlib.import_module(module_name)
            # Get a list of all the functions in the module
            function_list = [name for name in dir(module) if callable(getattr(module, name))]
            # Write the function names to the output file
            output_file.write(f"Functions in module {module_name_with_path}:\n")
            for function_name in function_list:
                output_file.write(f"{function_name}\n")
        # Check if the file is a folder
        elif os.path.isdir(file_path):
            # Recursively call this function on the subfolder
            list_functions_in_folder(file_path, output_file)

with open('FUNCTIONS.md', 'a') as f:
    list_functions_in_folder(".", f)