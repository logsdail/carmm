"""Crossover operation intended for bulk structures.
If you find this implementation useful in your work,
please consider citing:
    M. Van den Bossche, Henrik Gronbeck, B. Hammer,
    J. Chem. Theory Comput., doi:10.1021/acs.jctc.8b00039
in addition to the papers mentioned in the docstrings."""
import numpy as np
from random import random, randrange
from ase import Atoms
from ase.geometry import find_mic
from ase.build import niggli_reduce
from ase.ga.utilities import (atoms_too_close, atoms_too_close_two_sets,
                              gather_atoms_by_tag)
from ase.ga.offspring_creator import OffspringCreator


class Positions(object):
    """Helper object to simplify the pairing process.

    Parameters:

    scaled_positions: positions in scaled coordinates
    cop: center-of-positions (also in scaled coordinates)
    symbols: string with the atomic symbols
    distance: Signed distance to the cutting plane
    origin: Either 0 or 1 and determines which side of the plane
            the position should be at.
    """

    def __init__(self, scaled_positions, cop, symbols, distance, origin):
        self.scaled_positions = scaled_positions
        self.cop = cop
        self.symbols = symbols
        self.distance = distance
        self.origin = origin

    def to_use(self):
        """ Method which tells if this position is at the right side.
        """
        if self.distance > 0. and self.origin == 0:
            return True
        elif self.distance < 0. and self.origin == 1:
            return True
        else:
            return False


class CutAndSplicePairing(OffspringCreator):
    """ A cut and splice operator for bulk structures.

    For more information, see also:

    * `Glass, Oganov, Hansen, Comp. Phys. Comm. 175 (2006) 713-720`__

      __ https://doi.org/10.1016/j.cpc.2006.07.020

    * `Lonie, Zurek, Comp. Phys. Comm. 182 (2011) 372-387`__

      __ https://doi.org/10.1016/j.cpc.2010.07.048

    Parameters:

    blmin: dict
           The closest allowed interatomic distances on the form:
           {(Z, Z*): dist, ...}, where Z and Z* are atomic numbers.

    n_top: int or None
           The number of atoms to optimize (None = include all).

    p1: float or int between 0 and 1
        Probability that a parent is shifted over a random
        distance along the normal of the cutting plane.

    p2: float or int between 0 and 1
        Same as p1, but for shifting along the two directions
        in the cutting plane.

    minfrac: float or int between 0 and 1
             Minimal fraction of atoms a parent must contribute
             to the child.

    cellbounds: ase.ga.bulk_utilities.CellBounds instance
                Describing limits on the cell shape, see
                :class:`~ase.ga.bulk_utilities.CellBounds`.

    use_tags: boolean
              Whether to use the atomic tags to preserve
              molecular identity. Note: same-tag atoms are
              supposed to be grouped together.

    test_dist_to_slab: boolean
                       Whether also the distances to the slab
                       should be checked to satisfy the blmin.
    """

    def __init__(self, blmin, n_top=None, p1=1., p2=0.05, minfrac=None,
                 cellbounds=None, use_tags=False, test_dist_to_slab=True,
                 verbose=False):
        OffspringCreator.__init__(self, verbose)
        self.blmin = blmin
        self.n_top = n_top
        self.p1 = p1
        self.p2 = p2
        self.minfrac = minfrac
        self.cellbounds = cellbounds
        self.use_tags = use_tags
        self.test_dist_to_slab = test_dist_to_slab

        self.scaling_volume = None
        self.descriptor = 'CutAndSplicePairing'
        self.min_inputs = 1

    def update_scaling_volume(self, population, w_adapt=0.5, n_adapt=0):
        ''' Updates the scaling volume that is used in the pairing
        w_adapt: weight of the new vs the old scaling volume
        n_adapt: number of best candidates in the population that
                 are used to calculate the new scaling volume
        '''
        if not n_adapt:
            # take best 20% of the population
            n_adapt = int(round(0.2 * len(population)))
        v_new = np.mean([a.get_volume() for a in population[:n_adapt]])

        if not self.scaling_volume:
            self.scaling_volume = v_new
        else:
            volumes = [self.scaling_volume, v_new]
            weights = [1 - w_adapt, w_adapt]
            self.scaling_volume = np.average(volumes, weights=weights)

    def _get_pairing(self, a1, a2, direction=None, fraction=None):
        """
        Creates a child from two parents using the given cutting plane
        Does not check whether atoms are too close, but does return
        None if the generated structure lacks sufficient atoms
        from one of the parents (see self.minfrac).
        Assumes the parents have been 'topped' and checked for equal
        lengths, stoichiometries, and tags (if self.use_tags).
        direction: direction of the cutting surface normal (0, 1 or 2)
        fraction: fraction of the lattice vector along which
                  the cut is made.
        """
        N = len(a1)
        symbols = a1.get_chemical_symbols()
        pbc = a1.get_pbc()
        tags = a1.get_tags() if self.use_tags else np.arange(N)

        # Generate list of all atoms / atom groups:
        cell1 = a1.get_cell()
        cell2 = a2.get_cell()
        scalpos1 = a1.get_scaled_positions(wrap=False)
        scalpos2 = a2.get_scaled_positions(wrap=False)
        p1, p2, sym = [], [], []
        for i in np.unique(tags):
            indices = np.where(tags == i)[0]
            s = ''.join([symbols[j] for j in indices])
            sym.append(s)

            cop1 = np.mean(scalpos1[indices], axis=0)
            d1 = cop1[direction] - fraction
            p1.append(Positions(scalpos1[indices], cop1, s, d1, 0))

            cop2 = np.mean(scalpos2[indices], axis=0)
            d2 = cop2[direction] - fraction
            p2.append(Positions(scalpos2[indices], cop2, s, d2, 1))

        all_points = p1
        all_points.extend(p2)
        unique_sym = np.unique(sym)
        types = {s: sym.count(s) for s in unique_sym}

        # Sort these by chemical symbols:
        all_points.sort(key=lambda x: x.symbols, reverse=True)

        # For each atom type make the pairing
        unique_sym.sort()
        use_total = dict()
        for s in unique_sym:
            used = []
            not_used = []
            # The list is looked trough in
            # reverse order so atoms can be removed
            # from the list along the way.
            for i in reversed(range(len(all_points))):
                # If there are no more atoms of this type
                if all_points[i].symbols != s:
                    break
                # Check if the atom should be included
                if all_points[i].to_use():
                    used.append(all_points.pop(i))
                else:
                    not_used.append(all_points.pop(i))

            assert len(used) + len(not_used) == types[s] * 2

            # While we have too few of the given atom type
            while len(used) < types[s]:
                r = random()
                # origin = 0 => provides atoms if pos > fraction
                # origin = 1 => provides atoms if pos < fraction
                pick = 0 if r > fraction else 1
                interval = [0, fraction] if r < fraction else [fraction, 1]
                indices = []
                for index, point in enumerate(not_used):
                    cond1 = interval[0] <= point.cop[direction]
                    cond2 = point.cop[direction] <= interval[1]
                    if cond1 and cond2:
                        if point.origin != pick:
                            indices.append(index)
                if len(indices) == 0:
                    continue
                choice = randrange(0, len(indices))
                used.append(not_used.pop(choice))

            # While we have too many of the given atom type
            while len(used) > types[s]:
                # remove randomly:
                index = randrange(0, len(used))
                not_used.append(used.pop(index))

            use_total[s] = used

        n_tot = sum([len(ll) for ll in use_total.values()])
        assert n_tot == len(sym)

        # check if the generated structure contains atoms
        # from both parents:
        count1, count2 = 0, 0
        for x in use_total.values():
            count1 += sum([y.origin == 0 for y in x])
            count2 += sum([y.origin == 1 for y in x])

        nmin = 1 if self.minfrac is None else int(round(self.minfrac * N))
        if count1 < nmin or count2 < nmin:
            return None

        # pair the cells:
        if not self.scaling_volume:
            v_ref = 0.5 * (a1.get_volume() + a2.get_volume())
        else:
            v_ref = self.scaling_volume

        found = False
        while not found:
            r = random()
            newcell = r * cell1 + (1 - r) * cell2
            vol = abs(np.linalg.det(newcell))
            newcell *= (v_ref / vol)**(1. / 3)
            if self.cellbounds is not None:
                found = self.cellbounds.is_within_bounds(newcell)
            else:
                found = True

        # Construct the cartesian positions and reorder the atoms
        # to follow the original order
        newpos = []
        for s in sym:
            p = use_total[s].pop()
            c = cell1 if p.origin == 0 else cell2
            pos = np.dot(p.scaled_positions, c)
            cop = np.dot(p.cop, c)
            vectors, lengths = find_mic(pos - cop, c, pbc)
            newcop = np.dot(p.cop, newcell)
            pos = newcop + vectors
            for row in pos:
                newpos.append(row)

        newpos = np.reshape(newpos, (N, 3))
        num = a1.get_atomic_numbers()
        child = Atoms(numbers=num, positions=newpos, pbc=pbc, cell=newcell,
                      tags=tags)
        child.wrap()
        return child

    def get_new_individual(self, parents):
        """ The method called by the user that
        returns the paired structure. """
        f, m = parents

        indi = self.cross(f, m)
        desc = 'pairing: {0} {1}'.format(f.info['confid'],
                                         m.info['confid'])
        # It is ok for an operator to return None
        # It means that it could not make a legal offspring
        # within a reasonable amount of time
        if indi is None:
            return indi, desc
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid'],
                                        m.info['confid']]

        return self.finalize_individual(indi), desc

    def cross(self, a1, a2):
        """Crosses the two atoms objects and returns one"""

        if len(a1) != len(a2):
            raise ValueError('The two structures do not have the same length')

        N = len(a1) if self.n_top is None else self.n_top
        slab = a1[:len(a1) - N]
        a1 = a1[-N:]
        a2 = a2[-N:]

        if not np.array_equal(a1.numbers, a2.numbers):
            err = 'Trying to pair two structures with different stoichiometry'
            raise ValueError(err)

        if self.use_tags and not np.array_equal(a1.get_tags(), a2.get_tags()):
            err = 'Trying to pair two structures with different tags'
            raise ValueError(err)

        a1_copy = a1.copy()
        a2_copy = a2.copy()

        if self.cellbounds is not None:
            if not self.cellbounds.is_within_bounds(a1_copy.get_cell()):
                niggli_reduce(a1_copy)
            if not self.cellbounds.is_within_bounds(a2_copy.get_cell()):
                niggli_reduce(a2_copy)

        pos1_ref = a1_copy.get_positions()
        pos2_ref = a2_copy.get_positions()

        invalid = True
        counter = 0
        maxcount = 1000

        # Run until a valid pairing is made or 1000 pairings are tested.
        while invalid and counter < maxcount:
            counter += 1

            # Choose direction of cutting plane normal (0, 1, or 2):
            direction = randrange(3)

            # Randomly translate parent structures:
            for a, pos in zip([a1_copy, a2_copy], [pos1_ref, pos2_ref]):
                a.set_positions(pos)
                cell = a.get_cell()

                for i in range(3):
                    r = random()
                    cond1 = i == direction and r < self.p1
                    cond2 = i != direction and r < self.p2
                    if cond1 or cond2:
                        a.positions += random() * cell[i, :]

                if self.use_tags:
                    gather_atoms_by_tag(a)
                else:
                    a.wrap()

            # Perform the pairing:
            fraction = random()
            child = self._get_pairing(a1_copy, a2_copy, direction=direction,
                                      fraction=fraction)
            if child is None:
                continue

            # Verify whether the atoms are too close or not:
            invalid = atoms_too_close(child, self.blmin,
                                      use_tags=self.use_tags)
            if invalid:
                continue
            elif self.test_dist_to_slab:
                invalid = atoms_too_close_two_sets(slab, child, self.blmin)

        if counter == maxcount:
            return None

        return child
