class RotationBox():
    """

    Object intended to store functions for rotation

        Args:
            atoms_ads: Atoms object
                Contains positions of adsorbate.
            atoms_site: Atoms object
                Contains positions of the adsorption site.
            ads_idx: integer
                Index of the adsorbate atom bonding adsorbate and adsorption site.
            site_idx: integer
                Index of the site atom bonding adsorbate and adsorption site.
            bond_length: float
                Length of the desired bond between the site and adsorbate.
            neighb_idx: integer
                Index of neighbour used to form principle axes of rotation.
            lps: integer
                Number of lone pairs on the site atom.
            lp_idx: integer
                Index of the potential lone pair site (i.e., if two lone pair sites are present,
                controls which site the atom adsorbs to).
            cutoff_mult: float
                Constant factor multiplies species dependent cut-off radii used in
                carmm.analysis.neighbours to find the first nearest neighbour to the adsorption
                site (i.e., larger value potentially finds more neighbours).
        Returns:
            RotationBox: RotationBox Object
                Object used to adsorb a given adsorbate to a given site and rotate
                as desired using internal functions.

    """

    def __init__(self, atoms_ads, atoms_site, ads_idx, site_idx,
                 bond_length, neighb_idx=0, lps=1, lp_idx=1, cutoff_mult=1):

        self.atoms_ads = atoms_ads
        self.atoms_site = atoms_site

        self.ads_idx = ads_idx
        self.site_idx = site_idx
        self.bond_length = bond_length

        self.neighb_idx = neighb_idx
        self.lps = lps
        self.lp_idx = lp_idx

        self.cutoff_mult = cutoff_mult

        # Finds the normal (ie., the prospective bond to which the adsorbate is attached)
        self.site_norm = self.find_site_normal(self.atoms_site, self.site_idx, self.cutoff_mult)
        # Find the rotation axes the molecule hinges around.
        self.x_axis, self.y_axis, self.z_axis = self.find_adsorbate_rotation_axes()

    def rotate_and_place_adsorbate(self, rotation=[0, 0, 0]):
        """

        Place and rotates the adsorbate along a bond length defined by vector -
        vector calculated assuming VSEPR model with one remaining lone pair.
        Rotation takes places along the bonding direction. See documentation for
        rotation variable for directions of rotation.

        Args:
            atoms_ads: Atoms object
                Contains positions of adsorbate.
            atoms_site: Atoms object
                Contains positions of the adsorption site.
            bond_length: float
                Length of the bond for adsorption.
            ads_idx: integer
                Index of the adsorbate atom bonding adsorbate and adsorption site.
            site_idx: integer
                Index of the site atom bonding adsorbate and adsorption site.
            neighb_idx: integer
                Index of neighbour used to form principle axes of rotation.
            rotation: list
                Angles in degrees of rotation for absorbate.
                z - along ads and site bond, x - perpendicular to neighbour and bond, y - perpendicular to x and z.

        Returns:
            ads_and_site: numpy array, (:,3)
                Array of positions of the adsorbed species and site.
            atom_ads: numpy array, (:,3)
                Array of positions for the rotated adsorbed species only.

        """

        # Zeroes and rotates the molecule along the site normal
        junk, self.zeroed_adsorbate = self.place_adsorbate()

        self.zeroed_adsorbate.rotate(rotation[0], self.x_axis, center=self.zeroed_adsorbate.positions[self.ads_idx])
        self.zeroed_adsorbate.rotate(rotation[1], self.y_axis, center=self.zeroed_adsorbate.positions[self.ads_idx])
        self.zeroed_adsorbate.rotate(rotation[2], self.z_axis, center=self.zeroed_adsorbate.positions[self.ads_idx])

        self.ads_and_site = self.atoms_site + self.zeroed_adsorbate

    def rotate(self, rotation):

        import copy

        adsorbate = copy.deepcopy(self.zeroed_adsorbate)

        adsorbate.rotate(rotation[0], self.x_axis, center=adsorbate.positions[self.ads_idx])
        adsorbate.rotate(rotation[1], self.y_axis, center=adsorbate.positions[self.ads_idx])
        adsorbate.rotate(rotation[2], self.z_axis, center=adsorbate.positions[self.ads_idx])

        self.ads_and_site = self.atoms_site + adsorbate

        self.atoms_ads = copy.deepcopy(adsorbate)

    def place_adsorbate(self):
        """
        Place and rotates the adsorbate along a bond length defined by vector.
        Molecule rotated such that the center of positions for the molecule points along
        the bonding direction. Equivalent to rotate_and_place_adsorbate(...,rotation=[0,0,0]).

        Args:
            atoms_ads: Atoms object
                Contains positions of adsorbate.
            atoms_site: Atoms object
                Contains positions of the adsorption site.
            ads_idx: integer
                Index of the adsorbate atom bonding adsorbate and adsorption site.
            site_idx: integer
                Index of the site atom bonding adsorbate and adsorption site.
            bond_length: float
                Length of the bond for adsorption.

        Returns:
            ads_and_site: numpy array, (:,3)
                Array of positions of the adsorbed species and site.
            atom_ads: numpy array, (:,3)
                Array of positions for the rotated adsorbed species only.

        """

        import copy

        if len(self.atoms_ads.numbers) > 1:
            adsorbate_normal = self.find_adsorbate_normal(self.atoms_ads, self.ads_idx)

            # Returns the zeroed coordinates about atom[ads_idx].
            self.atoms_ads.positions = self.rotate_ads2site_vec(self.atoms_ads, self.site_norm,
                                                                adsorbate_normal)
        else:
            self.atoms_ads.positions = self.atoms_ads.positions - self.atoms_ads.positions

        # Places zeroed coordinates at the bonding position
        self.atoms_ads.positions = self.atoms_ads.positions + \
                                   (self.bond_length * self.site_norm) + self.atoms_site.positions[self.site_idx]

        self.ads_and_site = self.atoms_site + self.atoms_ads

        self.zeroed_adsorbate = copy.deepcopy(self.atoms_ads)

    def find_site_normal(self, atoms, index, cutoff_mult):
        """

        Returns a normalised vector specifying the direction of the adsorbate-adsorption
        site bond.

        Args:
            atoms: Atoms object
                Contains the atomic position of the site.
            index: integer
                The atomic index of the desired site atom.
            cutoff_mult: float
                Multiplier for the cutoff radii of neighbouring atoms (1 = natural cutoff radii).

        Returns:
            site_normal: numpy array, (3).
                Contains the normalised direction of the ads+site bond about (0,0,0).

        """
        from carmm.analyse.neighbours import neighbours
        from ase.geometry import find_mic
        from ase.neighborlist import natural_cutoffs
        import numpy as np

        assert self.lps < 3, "Lone pairs greater than 2 not yet implemented."
        assert self.lps != 0, "No valid adsorption site available."
        assert len(self.atoms_site) != 0, "Adsorbate site should have at least one other atom attached."

        cutoff = natural_cutoffs(atoms, cutoff_mult)
        neighbour_atoms, shell_list = neighbours(atoms, [index], 1, cutoff)

        vectors = atoms.positions[neighbour_atoms] - atoms.positions[index]

        if np.sum(atoms.get_cell().array) != 0.0:
            for v_idx, vector in enumerate(vectors):
                vectors[v_idx,:] = find_mic(vector, atoms.get_cell())[0]

        site_norm = np.sum(vectors, axis=0)

        site_norm = -site_norm / np.linalg.norm(site_norm)

        # Code and logic is sloppy and will be improved.
        if self.lps > 1:

            neighb_list, shell_list = neighbours(self.atoms_site, [self.site_idx], 1, cutoff)

            assert 1 < len(shell_list[1]) < 5, "Site either has too few or too many neighbours VSEPR."
            #        assert len(shell_list[1]) > 2, "Not implemented LPs above 2."

            if len(shell_list[1]) == 2:
                neighb1 = self.find_generic_normal(self.atoms_site,self.site_idx,shell_list[1][0])
                neighb2 = self.find_generic_normal(self.atoms_site,self.site_idx,shell_list[1][1])

                y_rot = np.cross(neighb1, neighb2)
                z_rot = np.cross(y_rot, site_norm)

                if self.lp_idx == 0:
                    theta = np.pi / 180 * -52
                elif self.lp_idx == 1:
                    theta = np.pi / 180 * 52

                rot_matrix = self.normal_rotation_matrix(theta, z_rot)

                site_norm = np.dot(rot_matrix, site_norm.T).T
                site_norm = site_norm / np.linalg.norm(site_norm)

        return site_norm

    def find_adsorbate_normal(self, atoms, index):
        """

        Returns a normalised vector specifying the direction of the desired bonding
        atom of the adsorbate, and the mean of the cartesian positions of the adsorbate
        molecule.

        Args:
            atoms: Atoms object
                Contains the atomic position of the adsorbate.
            index: integer
                The atomic index of the bonding atom.

        Returns:
            site_normal: numpy array, (3).
                Contains the normalised direction of bonding atom->center of positions
                about (0,0,0).

        """

        import numpy as np

        ads_com = atoms.get_positions().mean(axis=0)

        molecule_vector = (atoms.positions[index] - ads_com) / np.linalg.norm(atoms.positions[index] - ads_com)

        return molecule_vector

    def find_generic_normal(self, atoms, index1, index2):
        """
        Returns a normalised vector specifying the direction of a generic
        bond.

        Args:
            atoms: Atoms object
                Contains the atomic position of the adsorbate.
            index1: integer
                The atomic index of atom 1.
            index2: integer
                The atomic index of atom 2.

        Returns:
            site_normal: numpy array, (3).
                Contains the normalised direction of atom1 to atom2.
                about (0,0,0).

        """

        import numpy as np

        vector = atoms.positions[index2] - atoms.positions[index1]

        molecule_vector = (vector) / np.linalg.norm(vector)

        return molecule_vector

    def find_adsorbate_rotation_axes(self):
        """

        Find the x, y and z rotation axes, defined by the direction of the bond
        to the adsorbate (site_normal) and the direction of the bond to a given
        pre-existing neighbour atom. If only one around the adsorbate site is found,
        the next atom its neighbours are used.

        Args:
            atoms_site: Atoms object
                Contains positions of the adsorption site.
            site_idx: integer
                Index of the site atom bonding adsorbate and adsorption site.
            neighb_idx: integer
                Index of neighbour used to form principle axes of rotation.
            cutoff_mult: float
                Multiplier for the cutoff radii of neighbouring atoms (1 = natural cutoff radii).

        Returns:
            x_axis: numpy array, (3)
                x rotation axis (perpendicular to site bonds to the adsorbate and
                another neighbour atom).
            y_axis, numpy array, (3)
                y rotation axis (perpendicular given x and y rotation axes).
            z_axis, numpy array, (3)
                z rotation axis (defined by the site normal).

        """
        import numpy as np
        from carmm.analyse.neighbours import neighbours
        from ase.neighborlist import natural_cutoffs

        # Rotations performed on a RHS axis, with rotations about site normal z, x out plane wrt.
        # bond vector of another neighbour to the shell site and z, and y perpendicular to x and z.
        cutoff = natural_cutoffs(self.atoms_site, self.cutoff_mult)
        neighbour_atoms, shell_list = neighbours(self.atoms_site, [self.site_idx], 1, cutoff)

        # Clean up list for clarity.
        neighbour_atoms.remove(self.site_idx)
        assert len(neighbour_atoms) > 0, "Adsorption site does not have any neighbours. Cannot define rotation axes."

        # If the site only has one neighbour, the neighbour normal is unobtainable.
        # Next site is used to find an appropriate x and y rotation axis.
        if len(neighbour_atoms) == 1:
            new_site_idx = neighbour_atoms[0]
            neighbour_atoms, shell_list = neighbours(self.atoms_site, [new_site_idx], 1, cutoff)

            # Remove item of list with name of original, one neighbour atom.
            shell_list[1].remove(self.site_idx)

            # Order list for compatability reasons.
            shell_list[1].sort()

            vectors = self.atoms_site.positions[shell_list[1]] - self.atoms_site.positions[new_site_idx]
            neighb_vector = vectors[self.neighb_idx] / np.linalg.norm(vectors[self.neighb_idx])

        else:
            # Order list for compatability reasons.
            neighbour_atoms.sort()

            vectors = self.atoms_site.positions[neighbour_atoms] - self.atoms_site.positions[self.site_idx]
            neighb_vector = vectors[self.neighb_idx] / np.linalg.norm(vectors[self.neighb_idx])

        # Given the choice of two orthogonal vectors wrt the site normal is arbitrary, we choose the cross product
        # with respect to the vector from site to its neighbours, then the cross product with this newly generated
        # axis with the site normal.
        z_axis = self.site_norm
        x_axis = np.cross(neighb_vector, self.site_norm)
        y_axis = np.cross(x_axis, z_axis)

        return x_axis, y_axis, z_axis

    def rotate_ads2site_vec(self, atoms, site_vec, mol_vec):
        """

        Performs a rotation of the molecule onto the site normal.
        ROTATION MATRIX TAKEN FROM https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724

        Args:
            atoms: Atoms object
                Contains positions of adsorbate.
            atoms_idx: integer
                Index of the adsorbate atom bonding adsorbate and adsorption site.
            site_vec: numpy array, size=(3)
                vector defining the direction
            mol_vec: numpy array, size=(3)
                vector defining direction from the bonding atom to the mean of
                the atomic positions.

        Returns:
            rotated_positions: numpy array, size=(Nx3)
                Position of the atoms rotated so site_vec and mol_vec align.

        """

        import numpy as np

        theta = np.arccos(np.dot(mol_vec, site_vec))

        axis = np.cross(mol_vec, site_vec)

        axis = axis / np.linalg.norm(axis)

        rotation_matrix = -self.normal_rotation_matrix(theta, axis)

        zeroed_atomco = atoms.positions - atoms.positions[self.ads_idx]
        rotated_positions = np.dot(rotation_matrix, zeroed_atomco.T).T

        return rotated_positions

    def normal_rotation_matrix(self, theta, ax):
        """

        Returns a (3x3) rotation matrix for a rotation about a given
        axis by a theta value.

        Args:
            theta: float
                Angle of rotation in radians.
            ax: numpy array, size = (3)
                Direction of the rotation axis, centered on [0, 0, 0]

        Returns:
            rotation_matrix: numpy array, size = (3x3)
                Rotation matrix to be applied to atomic positions.

        """

        import numpy as np

        sina = np.sin(theta)
        cosa = np.cos(theta)

        onemincosa = 1.0 - cosa

        rotation_matrix = np.array([[(ax[0] * ax[0] * onemincosa) + cosa,
                                     (ax[1] * ax[0] * onemincosa) - (sina * ax[2]),
                                     (ax[2] * ax[0] * onemincosa) + (sina * ax[1])],
                                    [(ax[0] * ax[1] * onemincosa) + (sina * ax[2]),
                                     (ax[1] * ax[1] * onemincosa) + cosa,
                                     (ax[2] * ax[1] * onemincosa) - (sina * ax[0])],
                                    [(ax[0] * ax[2] * onemincosa) - (sina * ax[1]),
                                     (ax[1] * ax[2] * onemincosa) + (sina * ax[0]),
                                     (ax[2] * ax[2] * onemincosa) + cosa]])

        return rotation_matrix
