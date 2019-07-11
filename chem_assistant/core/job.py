from .atom import Atom
from .molecule import Molecule
from .settings import *
from .sc import Supercomp
from .utils import (sort_elements, write_xyz)

from os.path import (basename, dirname, join, exists)
from os import(mkdir, chdir, getcwd, system, walk, listdir)
import sys

__all__ = ['Job']

class Job:
    """Base class for any input file for a computational chemistry calculation- ab initio or
molecular dynamics. This class also creates job files in the same directory as the class is called. 

    Instances of this class have the following attributes:
    * ``using`` -- coordinates of chemical system, in xyz format

    """
    def __init__(self, using = None, run_dir = None, frags_in_subdir = False, **kwargs): # allows for fmo=True, even if nothing done with the arguments
        if using is not None:
            self.molecule_name = using
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
                          'mon': 'mon',
                          'monarch': 'mon',
                          'mas': 'mas',
                          'massive': 'mas',
                          'm3': 'mas',
                          'gaia': 'gaia',
                          'stm': 'stm',
                          'stampede': 'stm'}
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

    def write_file(self, data, filetype):
        """Writes the generated PSI4 input/jobs to a file. If no filename is passed when the class is instantiated, the name of the file defaults to the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (freq). 

        NOTE: Must pass data as a string, not a list!""" 
        with open(f"{self.base_name}.{filetype}", "w") as f:
            f.write(data)

    def get_job_template(self):
        job_file = self.find_job()
        with open(job_file) as f:
            job = f.read()       
            return job

    def create_job(self):
        """Returns the relevant job template as a list, then performs the necessary modifications. After, the job file is printed in the appropriate directory."""
        jobfile = self.get_job_template()
        # modify
        if str(self.sc) == 'mgs':
            jobfile = self.change_mgs_job(jobfile)
            jobfile = jobfile.replace('name', f'{self.base_name}') 
        elif str(self.sc) == 'rjn':
            jobfile = self.change_rjn_job(jobfile)
            jobfile = jobfile.replace('name', f'{self.base_name}') 
        elif str(self.sc) == 'mas':
            jobfile = jobfile.replace('base_name', f'{self.base_name}') 
        elif str(self.sc) == 'mon':
            jobfile = jobfile.replace('base_name', f'{self.base_name}') 
        elif str(self.sc) == 'stm':
            jobfile = jobfile.replace('base_name', f'{self.base_name}') 
        self.write_file(jobfile, filetype='job')
