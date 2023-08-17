import os
from fnmatch import fnmatch
from ase.io import read

class CalculationHelper:
    def __init__(self,
                 calc_type,
                 parent_dir,
                 filename,
                 restart=True,
                 verbose=True,
                 prev_calcs=None,
                 interpolation=None):

        self.calc_type = calc_type
        self.parent_dir = parent_dir
        self.filename = filename
        self.restart = restart
        self.verbose = verbose
        self.counter = 0
        self.prev_calcs = None
        self.interpolation = None

    def _find_restart(self):
        restart_found = False
        initial = None
        while not restart_found and self.counter > 0:
            subdirectory_name_previous = f"{self.calc_type}_{self.filename}_{self.counter - 1}"

            # Handle different types of calculations
            if self.calc_type == "Opt":
                opt_restarts = 0
                while os.path.exists(os.path.join(subdirectory_name_previous, f"{self.counter - 1}_{self.filename}_{opt_restarts}.traj")):
                    opt_restarts += 1

                traj_name = os.path.join(subdirectory_name_previous, f"{self.counter - 1}_{self.filename}_{opt_restarts - 1}.traj")

                while os.path.exists(traj_name):
                    if os.path.getsize(traj_name):
                        if self.verbose:
                            print(f"Restarting calculation from {traj_name}")
                        initial = read(traj_name)
                        restart_found = True
                        break

                    traj_name = os.path.join(subdirectory_name_previous, f"{self.counter - 1}_{self.filename}_{opt_restarts - 1}.traj")
                    opt_restarts -= 1

            self.counter -= 1

        return initial

    def restart_setup(self):
        supported_calc_types = ["Opt", "Vib", "TS", "Charges"]
        assert self.calc_type in supported_calc_types

        while f"{self.calc_type}_{self.filename}_{self.counter}" in [fn for fn in os.listdir(self.parent_dir) if fnmatch(fn, f"{self.calc_type}*")]:
            self.counter += 1

        if self.counter > 0 and self.restart:
            initial = self._find_restart()
            self.counter += 1 # important to ensure new calculation begins in a new folder
        else:
            initial = None

        subdirectory_name = f"{self.calc_type}_{self.filename}_{self.counter}"
        out = f"{self.counter}_{self.filename}.out"



        return self.counter, out, subdirectory_name, initial
