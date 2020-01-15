from .atom import Atom
from .molecule import Molecule
from .settings import *
from .sc import Supercomp
from .utils import (sort_elements, write_xyz)

from os.path import (basename, dirname, join, exists)
from os import (mkdir, chdir, getcwd, system, walk, listdir)
import sys

__all__ = ['Job']

class Job:
    """Base class for any input file for a computational chemistry calculation- ab initio or
molecular dynamics. This class also creates job files in the same directory as the class is called. 

    Instances of this class have the following attributes:
    * ``using`` -- coordinates of chemical system, in xyz format

    """

    SLURM_HOSTS = ('stm', 'mas', 'mon')
    PBS_HOSTS = ('rjn', 'gadi')

    def __init__(self, using = None, run_dir = None, frags_in_subdir = False, user_settings=None,
    bonds_to_split = None, **kwargs): 
        # allows for fmo=True, even if nothing done with the arguments
        # pass on grouping/splitting to the base Molecule class
        if using is not None:
            self.molecule_name = using
            if user_settings is not None and 'grouped' in user_settings.keys():
                self.mol = Molecule(using, group = user_settings.grouped, bonds_to_split=bonds_to_split)
            else:
                self.mol = Molecule(using, bonds_to_split=bonds_to_split)

    def __repr__(self):
        return f'{self.__class__.__name__}: {self.mol.xyz}'

    __str__ = __repr__

    def get_sc(self):
        """
        Attempts to read in the supercomputer from a user defined `sett.supercomp`,
        or uses the hostname of the cluster to decide otherwise
        """
        if hasattr(self, 'merged'):
            meta = self.merged
        else:
            meta = self.defaults
        if 'supercomp' in meta.keys(): #user has to define a supercomp
            user_sc = meta.supercomp
            supercomps = {'rjn': 'rjn',
                          'raijin': 'rjn',
                          'gadi': 'gadi',
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
            self.sc = Supercomp().supercomp 

    def find_job(self, dft=False):
        """
        Returns the relevant job template. If a GAMESS job is for a dft 
        calculation, gamess_{self.sc}_dft.job will be called"""
        self.get_sc()
        package = sys.modules[self.__class__.__module__].__file__.split('/')[-1][:-3]
        # Returns gamess from GamessJob, psi from PsiJob etc...
        if dft:
            job = f"{package}_{self.sc}_dft.job"
        else:
            job = f"{package}_{self.sc}.job"
        dir_name = dirname(__file__) 
        templates = join(dir_name, '..', 'templates')
        job_file = join(templates, job)
        return job_file

    def find_charge_and_mult(self):
        """
        Changes charge and multiplicity unless user defines values in a
        settings file. In that case, the user-defined charge and 
        multiplicity are used.
        """
        user_assigned_charge=False
        user_assigned_mult=False
        if hasattr(self, 'user_settings') and self.user_settings != {}:
            if 'charge' in self.user_settings['input'].keys():
                user_assigned_charge=True
            if 'mult' in self.user_settings['input'].keys():
                user_assigned_mult=True
        if not user_assigned_charge:
            self.input.charge = self.mol.overall_charge
        if not user_assigned_mult:
            self.input.mult = self.mol.overall_mult        

    def write_file(self, data, filetype):
        """Writes the generated input/jobs to a file. If no filename is passed when the class is instantiated, the name of the file defaults to the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (freq). 

        NOTE: Must pass data as a string, not a list!""" 
        with open(f"{self.base_name}.{filetype}", "w") as f:
            f.write(data)

    def get_job_template(self, dft=False):
        job_file = self.find_job(dft=dft)
        with open(job_file) as f:
            job = f.read()       
            return job

    # def create_job(self):
    #     """Returns the relevant job template as a list, then performs the necessary modifications. After, the job file is printed in the appropriate directory."""
    #     jobfile = self.get_job_template()
    #     # modify
    #     if str(self.sc) == 'mgs':
    #         jobfile = self.change_mgs_job(jobfile)
    #         jobfile = jobfile.replace('name', f'{self.base_name}') 
    #     elif str(self.sc) == 'rjn':
    #         jobfile = self.change_rjn_job(jobfile)
    #         jobfile = jobfile.replace('name', f'{self.base_name}') 
    #     elif str(self.sc) == 'mas':
    #         jobfile = jobfile.replace('base_name', f'{self.base_name}') 
    #     elif str(self.sc) == 'mon':
    #         jobfile = jobfile.replace('base_name', f'{self.base_name}') 
    #     elif str(self.sc) == 'stm':
    #         jobfile = jobfile.replace('base_name', f'{self.base_name}') 
    #     self.write_file(jobfile, filetype='job')
