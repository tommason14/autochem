from .periodic_table import PeriodicTable as PT
import math
import numpy as np

__all__ = ["Atom"]


class Atom:
    """A class representing a atom in 3 dimensional euclidean space.
    An instance has the following attributes:

    * ``atnum`` -- atomic symbol, equal to zero for a dummy atom
    * ``coords`` -- tuple of x,y,z coordinates
    * ``bonds`` -- list of bonds that this atom is a part of
    * ``mol`` -- molecule this atom is a part of. Assigned programmatically when a molecule is separated using the *mol.separate* method, or can be assigned manually if building up a molecule from scratch

    Access these properties directly:
    * ``x``, ``y``, ``z`` -- for atom coordinates
    * ``symbol`` -- read or write the atom symbol directly

    >>> a = Atom('H', coords = (1,2,3))

    """

    def __init__(self, symbol=None, atnum=0, coords=None, mol=None, bonds=None):
        if symbol is not None:
            self.symbol = symbol
            self.atnum = PT.get_atnum(self)
        else:
            self.atnum = atnum
            self.symbol = PT.get_symbol(self)

        self.mass = PT.get_mass(self)
        self.mol = mol
        self.bonds = bonds or []
        self.connected_atoms = []
        self.h_bonded_to = []
        self.fragment = None

        if coords is None:
            self.coords = [0, 0, 0]
        elif len(coords) == 3:
            self.coords = [float(i) for i in coords]
        else:
            raise TypeError("Atom: Invalid coordinates given")

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]

    @property
    def z(self):
        return self.coords[2]

    def __repr__(self):
        """Unambiguous representation of an |Atom| instance"""
        if self.fragment is not None:
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Index: {self.index} Mol: {self.fragment}"
        if hasattr(self, "index") and not hasattr(self, "mol"):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Index: {self.index} "
        if hasattr(self, "index") and hasattr(self, "mol"):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Index: {self.index} Mol: {self.mol}"
        elif hasattr(self, "number"):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Mol: {self.mol} Atom: {self.number}"
        elif hasattr(self, "index") and len(self.h_bonded_to) > 0:
            h_bonded = [
                (atom.symbol, {"mol": atom.mol, "atom": atom.index})
                for atom in self.h_bonded_to
            ]
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Mol: {self.mol} Index: {self.index} Number: {self.number} H-Bonds: {h_bonded}"
        return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}"

    def __iter__(self):
        """Iterates through coordinates when called"""
        return iter(self.coords)

    def translate(self, vector):
        """Move atom in space by passing a vector in angstroms"""
        self.coords = [i + j for i, j in zip(self, vector)]

    def move_to(self, vector):
        """Move atom in space to the values, in angstroms, given in this vector. The vector passed represents a point in euclidean space"""
        self.coords = [i for i in vector]

    def distance_to(self, vector):
        """Measure the distance between the atom and a point in space, given as a vector in angstroms"""
        # pythagoras in 3D
        dist = 0.0
        for i, j in zip(self, vector):
            dist += (i - j) ** 2
        return dist ** 0.5

    def vector_to(self, point):
        """Returns a vector from the atom to a given point, in angstroms"""
        return [(i - j) for i, j in zip(point, self)]

    def angle_between(self, pos1, pos2):
        """Returns an angle between positions 1 and 2 in degrees, with this atom lying at the centre"""
        # dot product, angle = cos^-1([vec(a).vec(b)] / [dist(a) * dist(b)])
        num = np.dot(self.vector_to(pos1), self.vector_to(pos2))
        denom = self.distance_to(pos1) * self.distance_to(pos2)
        return math.acos(num / denom) * (180 / math.pi)

    def as_xyz(self, dps=5, end_of_line="\n"):
        """
        Return atom in xyz format: symbol x y z. Can also give an optional 
        end of line character such as a space, as well as specify how many 
        decimal places the coordinates should be given to.
        """
        return f"{self.symbol:5s} {self.x:>13.{dps}f} {self.y:>13.{dps}f} {self.z:>13.{dps}f}{end_of_line}"

    __str__ = __repr__
