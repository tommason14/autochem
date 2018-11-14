from .atom import Atom
from .molecule import Molecule
from .settings import *

from os.path import (basename, dirname, join)
import sys

__all__ = ['Job']

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

    # get_sc() and create_job() are meaningless when not used in either interfaces/gamess.py or
    # interfaces/psi.py

    def get_sc(self):
        if hasattr(self, 'merged'):
            meta = self.merged
        else:
            meta = self.defaults
        if 'supercomp' in meta.keys(): #user has to define a supercomp
            user_sc = meta.supercomp
            supercomps = {'rjn': 'rjn',
                          'raijin': 'rjn',
                          'mgs': 'mgs',
                          'magnus': 'mgs',
                          'gaia': 'gaia'}
            try:
                self.sc = supercomps[user_sc]
            except:
                raise AttributeError('Please enter a different, more specific string for the supercomputer- or remove the declaration and let the program decide.')
        else:
            self.sc = Supercomp()

    def find_job(self):
        """Returns the relevant job template"""
        self.get_sc()
        package = sys.modules[self.__class__.__module__].__file__.split('/')[-1][:-3]
        # package is the basename of the file containing the class GamessJob or PsiJob, which
        # inherit from Job- so we're only using this in inherited classes - returns gamess or psi
        job = f"{package}_{self.sc}.job"
        dir_name = dirname(__file__) #  (dir of job.py = core)
        templates = join(dir_name, '..', 'templates')
        job_file = join(templates, job)
        return job_file
