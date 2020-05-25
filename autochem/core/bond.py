__all__ = ["Bond"]


class Bond:
    """A class representing the bond between two atoms

    An instance of the class has the following attributes:

    * ``atom1`` and ``atom2`` -- two instances of |Atom| linked by this bond
    * ``order`` -- bond order; integer or floating point
    * ``mol`` -- |Molecule| this bond belongs to

    Newly created bonds are not added to ``atom1.bonds`` or ``atom2.bonds``- this information is added via the ``Molecule.add_bond`` method. """

    AROMATIC_ORDER = 1.5

    def __init__(self, atom1=None, atom2=None, order=1, mol=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.mol = mol

    def __repr__(self):
        """Returns a string representation of the bond"""
        return (
            f"({str(self.atom1).strip()})--{self.order:1.1f}--({str(self.atom2).strip()})"
        )

    __str__ = __repr__

    def __iter__(self):
        """Iterate over the bonded atoms, ``atom1`` first"""
        yield self.atom1
        yield self.atom2

    def is_aromatic(self):
        """Check if bond is aromatic, by looking at bond order"""
        return self.order == AROMATIC_ORDER

    def length(self):
        """Returns bond length in angstroms"""
        return self.atom1.distance_to(self.atom2)

    def other_end(self, atom):
        """Returns the atom on the other end of the bond"""
        if atom is self.atom1:
            return self.atom2
        elif atom is self.atom2:
            return self.atom1
        else:
            return AttributeError("Bond.other_end: invalid atom used")

    def resize(self, atom, length):
        """Change length of bond to length given, in angstroms"""
        moving_atom = self.other_end(atom)
        moving_atom.translate(tuple(i for i in moving_atom.vector_to(atom)))
