# Author: Igor Kowalec

class CalculationHelper:
    """
    A class for aiding setup of DFT calculations in an ASE/i-Pi/FHI-aims configuration using carmm.run.workflows.react.
    """

    def __init__(self,
                 calc_type: str,
                 parent_dir: str,
                 filename: str,
                 restart=True,
                 verbose=True):
        '''
        Args:
            calc_type: str
                Define calculation type to choose the setup and restart routine. Valid choices are currently:
                "Opt", "Charges", "TS". Note that vibration calculations are not currently governed by the
                Calculation Helper.
            parent_dir: str
                The parent directory of each CalculationHelper instance is be controlled be setting the parent_dir
                excplicitly, this is subsequently used for detecting previous calculations in restart_setup().
            filename: str
                String containing the name of the file to be used for detecting restarts and setting up working
                directories.
            restart: bool
                Attempt to restart the calculation by looking through detected previous working directories and
                reading the relevant files. If False, a new calculation will be restarted from scratch.
            verbose: bool
                Control the verbosity of the calculation setup by turning on (True) or off (False).
        '''

        self.calc_type = calc_type
        self.parent_dir = parent_dir
        self.filename = filename
        self.restart = restart
        self.verbose = verbose
        self.counter = 0

        '''Placeholders for properties retrieved from TS calculations'''
        self.prev_calcs = None
        self.interpolation = None

    def _find_restart(self):
        '''
        A function searching for previous converged calculations based on the calculation type in self.calc_type and
        the naming convention based on self.filename.

        Returns:
            For self.calc_type in ["Opt", "Charges"] returns an Atoms object
            For self.calc_type == "TS" returns a list of Atoms object
        '''
        import os
        from ase.io import read

        restart_found = False
        initial = None
        while not restart_found and self.counter > 0:
            subdirectory_name_previous = f"{self.calc_type}_{self.filename}_{self.counter - 1}"

            # Handle different types of calculations
            # Note that Vib does not need to be handled here as restarts are efficient in ase.vibrations
            if self.calc_type == "Opt":
                opt_restarts = 0
                while os.path.exists(os.path.join(subdirectory_name_previous,
                                                  f"{self.counter - 1}_{self.filename}_{opt_restarts}.traj")):
                    opt_restarts += 1

                traj_name = os.path.join(subdirectory_name_previous,
                                         f"{self.counter - 1}_{self.filename}_{opt_restarts - 1}.traj")

                while os.path.exists(traj_name) and not restart_found:
                    if os.path.getsize(traj_name):
                        if self.verbose:
                            print(f"Restarting calculation from {traj_name}")
                        initial = read(traj_name)
                        restart_found = True
                    if not restart_found:
                        traj_name = os.path.join(subdirectory_name_previous,
                                                 f"{self.counter - 1}_{self.filename}_{opt_restarts - 1}.traj")
                        opt_restarts -= 1

            elif self.calc_type == "Charges":
                traj_name = os.path.join(subdirectory_name_previous,
                                         f"{self.filename}_{self.calc_type.lower()}.traj")
                if os.path.exists(traj_name):
                    if os.path.getsize(traj_name):
                        if self.verbose:
                            print(f"Restarting calculation from {traj_name}")
                        initial = read(traj_name)
                        restart_found = True

            elif self.calc_type == "TS":

                initial = [None, None]
                traj_name = f"{subdirectory_name_previous}/ML-NEB.traj"
                if os.path.exists(traj_name):
                    print("TS search already converged at", traj_name)
                    minimum_energy_path = read(f"{traj_name}@:")
                    initial[0] = minimum_energy_path
                    restart_found = True
                else:
                    minimum_energy_path = read(f"{subdirectory_name_previous}/last_predicted_path.traj@:")
                    initial[1] = minimum_energy_path
                    restart_found = True

            # Point to the last folder
            self.counter -= 1

        # Important to ensure new calculation begins in a new folder - adjust outside the while loop
        self.counter += 1

        return initial

    def restart_setup(self):
        """
        The function for initializing the calculation based on self.calc_type and requested self.restart.
        It returns the counter required for setting up the working directory, the FHI-aims output filename
        following naming convention, the name of the working subdirectory and the previously converged structure(s).

        Returns:
            self.counter: int, out: str, subdirectory_name: str, initial: Atoms or list of Atoms objects

        """
        from fnmatch import fnmatch
        import os

        supported_calc_types = ["Opt", "Vib", "TS", "Charges"]
        assert self.calc_type in supported_calc_types

        while f"{self.calc_type}_{self.filename}_{self.counter}" in [
                fn for fn in os.listdir(self.parent_dir) if fnmatch(fn, f"{self.calc_type}*")]:
            self.counter += 1

        if self.counter > 0 and self.restart:
            initial = self._find_restart()
        else:
            initial = None

        subdirectory_name = f"{self.calc_type}_{self.filename}_{self.counter}"
        out = f"{self.counter}_{self.filename}.out"

        return self.counter, out, subdirectory_name, initial
