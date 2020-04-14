# TODO: Duplicate this for spin.
def extract_mulliken_charge(fn, natoms):
    '''
    Function to extract and return the Mulliken charges from an FHI-aims output.
    A list of relative charges is returned, with +q meaning charg depletion and -q meaning charge accummulation

    Parameters:

    fn: string
        Filename from which the Mulliken data should be extracted
    natoms: int
        Number of atoms in the calculation

    '''

    with open(fn, 'r') as f:
        output = f.readlines()

    # Focus on just the Mulliken data
    mulliken_line = 0
    output_line = len(output) - 1
    while output_line and not mulliken_line:
        if "Starting Mulliken Analysis" in output[output_line]:
            mulliken_line = output_line
            # 8 lines from text label to start of data (spin-paired)
            mulliken_line += 8
        else:
            output_line -= 1

    mulliken_data = [q.split()[3] for q in output[mulliken_line:mulliken_line + natoms]]

    return mulliken_data

def parse_mulliken_file(lines):
    '''

    :param lines:
    :return:
    '''

    natoms = 1
    nspin = 1
    nkpts = 1
    nstates = 1
    # Prescan Mulliken file for natoms, nspin, nkpts so we can then digest data.
    for line in lines:
        words = line.split()
        if (len(words)) > 0:
            if words[0] == "Atom":
                # Outer loop for Mulliken is over each atom
                current_atom = int(words[2].replace(":", ""))
                if current_atom > natoms:
                    natoms = current_atom
            elif words[0] == "Spin":
                # Middle loop is over spins. Keyword only present if spin is present
                nspin = 2
            elif words[0] == "k":
                # Inside loop is over k-points.
                current_k_point = int(words[3].replace(":",""))
                if current_k_point > nkpts:
                    nkpts = current_k_point
            else:
                try:
                    current_state = int(words[0])
                    if current_state > nstates:
                        nstates = current_state
                except ValueError:
                    pass

    # Setup all data structures containing all information for each atom
    md = MullikenData(natoms, nspin, nkpts, nstates)

    current_atom = 1
    current_spin = 0
    current_k_point = 1

    # Spin keywords are different depending on whether included in the calculation
    # For spin-paired systems, we won't hit a keyword and so need to force increment this value
    if nspin == 1:
        current_spin = 1

    for line in lines:
        words = line.split()
        if len(words) > 0:
            if words[0] == "Atom":
                # Outer loop for Mulliken is over each atom
                current_atom = int(words[2].replace(":", ""))
            elif words[0] == "Spin":
                current_spin += 1
                if current_spin > 2:
                    current_spin = 1
            elif words[0] == "k":
                # Inside loop is over k-points.
                current_k_point = int(words[3].replace(":",""))
                if len(words) > 9:
                    md.atoms[current_atom-1].spin[current_spin-1].kpts[current_k_point-1].weight = float(words[10])
            else:
                try:
                    current_state = int(words[0])
                    md.atoms[current_atom-1].spin[current_spin-1].kpts[current_k_point-1].energies[
                        current_state-1] =  float(words[1])
                    md.atoms[current_atom-1].spin[current_spin-1].kpts[current_k_point-1].occupancies[
                        current_state-1] = float(words[2])
                    md.atoms[current_atom-1].spin[current_spin-1].kpts[current_k_point-1].all_mulliken[
                        current_state-1] = float(words[3])
                    for i in range(len(words) - 4):
                        md.atoms[current_atom-1].spin[current_spin-1].kpts[current_k_point-1].orbitals[i][current_state-1] = float(words[i+4])
                except ValueError:
                    pass

    return md

class Kpt:
    def __init__(self, nstates):
        self.weight = 1.0
        self.energies = [ 0.0 for i in range(nstates) ]
        self.occupancies = [ 0.0 for i in range(nstates) ]
        self.all_mulliken = [ 0.0 for i in range(nstates) ]
        self.orbitals = [ [ 0.0 for i in range(nstates) ] for j in range(5) ] # Up to l = 4, which is g.

class Spin:
    def __init__(self, nkpts, nstates):
        self.kpts = [ Kpt(nstates) for i in range(nkpts) ]

class Atom:
    def __init__(self, nspin, nkpts, nstates):
        self.spin = [ Spin(nkpts, nstates) for i in range(nspin) ]

class MullikenData:
    def __init__(self, natoms, nspin, nkpts, nstates):
        self.atoms = [ Atom(nspin, nkpts, nstates) for i in range(natoms) ]
        self.homo = None

    def get_natoms(self):
        return len(self.atoms)

    def get_nspin(self):
        return len(self.atoms[0].spin)

    def get_nkpts(self):
        return len(self.atoms[0].spin[0].kpts)

    def get_nstates(self):
        return len(self.atoms[0].spin[0].kpts[0].energies)

    def get_homo(self):
        if self.homo is None:
            # Arbitrarily set the HOMO to the first state we can sample
            self.homo = self.atoms[0].spin[0].kpts[0].energies[0]
            # Iterate through all atoms, spins and k-points to get the true HOMO
            for atom in range(self.get_natoms()):
                for sp in range(self.get_nspin()):
                    for kpt in range(self.get_nkpts()):
                        # Get first value for this k-point where the occupancies are below 0.5
                        # We then subtract one from this to get the last "occupied" state
                        # Note: This doesn't deal well with delocalised charge over degenerate states.
                        # TODO: Configure so it'll test subsequent states for degeneracy and split occupancy
                        res = next(x for x, val in enumerate(self.atoms[atom].spin[sp].kpts[kpt].occupancies)
                                if val < 0.5) - 1
                        if self.atoms[atom].spin[sp].kpts[kpt].energies[res] > self.homo:
                            self.homo = self.atoms[atom].spin[sp].kpts[kpt].energies[res]

        # Even though we calculate the HOMO, for periodic systems the x-axis is shifted so
        # the HOMO is zero. Therefore return zero to match the x-axis returned.
        if self.get_nkpts() > 1:
            return 0.0
        else:
            return self.homo

    def get_all_plot_data(self):
        return self.get_plot_data(atoms=range(self.get_natoms()),
                                  spin=range(self.get_nspin()),
                                  kpts=range(self.get_nkpts()),
                                  angular='all')

    def get_s_plot_data(self, atoms=None, spin=None, kpts=None):
        if atoms is None:
            atoms = range(self.get_natoms())
        if spin is None:
            spin = range(self.get_nspin())
        if kpts is None:
            kpts = range(self.get_nkpts())
        return self.get_plot_data(atoms, spin, kpts, 's')

    def get_p_plot_data(self, atoms=None, spin=None, kpts=None):
        if atoms is None:
            atoms = range(self.get_natoms())
        if spin is None:
            spin = range(self.get_nspin())
        if kpts is None:
            kpts = range(self.get_nkpts())
        return self.get_plot_data(atoms, spin, kpts, 'p')

    def get_d_plot_data(self, atoms=None, spin=None, kpts=None):
        if atoms is None:
            atoms = range(self.get_natoms())
        if spin is None:
            spin = range(self.get_nspin())
        if kpts is None:
            kpts = range(self.get_nkpts())
        return self.get_plot_data(atoms, spin, kpts, 'd')

    def get_f_plot_data(self, atoms=None, spin=None, kpts=None):
        if atoms is None:
            atoms = range(self.get_natoms())
        if spin is None:
            spin = range(self.get_nspin())
        if kpts is None:
            kpts = range(self.get_nkpts())
        return self.get_plot_data(atoms, spin, kpts, 'f')

    def get_plot_data(self, atoms, spin, kpts, angular, ymin=-20, ymax=+20, npoints=1000, variance=0.02):

        import numpy as np
        from scipy.stats import norm

        # Dictionary to simplify pulling out the angular decomposition
        angular_momenta = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

        # Make sure ymin and ymax are sane
        if ymin > ymax:
            temp = ymax
            ymax = ymin
            ymin = temp
        # Plotting variables
        data = [ [ 0.0 ] * npoints ] * (len(spin))
        x = np.linspace(ymin, ymax, npoints)
        sigma = np.sqrt(variance)

        # Check whether we have the HOMO; if not, calculate.
        if self.homo is None:
            homo = self.get_homo()

        # Collect information in data object. Note two sets of results, for spin up and down.
        for atom in atoms:
            for sp in spin:
                for kpt in kpts:
                    for e in range(self.get_nstates()):
                        energy = self.atoms[atom].spin[sp].kpts[kpt].energies[e]

                        if energy > ymin and energy < ymax:
                            if angular == 'all':
                                data[sp] += self.atoms[atom].spin[sp].kpts[kpt].all_mulliken[e] * \
                                            self.atoms[atom].spin[sp].kpts[kpt].weight * \
                                            norm.pdf(x, energy, sigma)
                            else:
                                for n in angular:
                                    data[sp] += self.atoms[atom].spin[sp].kpts[kpt].orbitals[angular_momenta[n]][e] * \
                                                self.atoms[atom].spin[sp].kpts[kpt].weight * \
                                                norm.pdf(x, energy, sigma)

        if self.get_nkpts() > 1:
            return x-self.homo, data
        else:
            return x, data

    def get_data_integrity(self):
        x, md = self.get_all_plot_data()
        x, md_spdfg = self.get_plot_data(range(self.get_natoms()),
                                         range(self.get_nspin()),
                                         range(self.get_nkpts()),'spdfg')

        total_diff = 0.0
        total_md = 0.0
        for sp in range(self.get_nspin()):
            #print(abs(sum(md[sp] - md_spdfg[sp])), "( " + str(100*abs(sum(md[sp] - md_spdfg[sp]))/sum(md[sp])) + "% )" )
            total_diff += abs(sum(md[sp] - md_spdfg[sp]))
            total_md += sum(md[sp])

        # Checks that the difference of integrals of data are less than 0.001 %
        try:
            assert(total_diff/total_md < 1e-5)
            return True
        except:
            print("Data integrity seems compromised - check your data rigorously before further use")
            return False

    def get_graph_xlabel(self):
        if self.get_nkpts() > 1:
            xlabel = '$\epsilon - \epsilon_{f}$ (eV)'
        else:
            xlabel = '$\epsilon$ (eV)'
        return xlabel
