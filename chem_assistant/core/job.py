from .atom import Atom
from .molecule import Molecule
from .settings import *

class Job:
    """Base class for any input file for a computational chemistry calculation- ab initio or
molecular dynamics. This class also creates job files in the same directory as the class is called. 

    Instances of this class have the following attributes:
    * ``using`` -- coordinates of chemical system, in xyz format

    """
    def __init__(self, using = None):
        if using is not None:
            self.mol = Molecule(using)
            #self.mol.separate() # in the "user side" file
            #self.mol.gamess_format() # in the "user side" file
            # for frag in self.mol.fragments: Input_type(using=f'fragments/{frag}') -> no fmo though


