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
        if run_dir is not None:
            self.made_run_dir = True
        else:
            self.made_run_dir = False

    # get_sc() and create_job() are meaningless when not used in either interfaces/gamess.py or
    # interfaces/psi.py

    def title(self):
        if '/' in self.molecule_name:
            self.title = self.molecule_name.split('/')[-1][:-4] #say using = ../xyz_files/file.xyz --> 
        else:
            self.title = self.molecule_name[:-4]

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

    def write_file(self, data, filetype):
        """Writes the generated PSI4 input/jobs to a file. If no filename is passed when the class is instantiated, the name of the file defaults to the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (freq). 

        NOTE: Must pass data as a string, not a list!""" 
        with open(f"{self.base_name}.{filetype}", "w") as f:
            f.write(data)


    def place_files_in_dir(self):
            """Move input and job files into a directory named with the input name (``base_name``) i.e.
            moves opt.inp and opt.job into a directory called ``opt``."""
            if self.is_complex:
                if not exists(join(self.base_name, 'complex')):
                    mkdir(join(self.base_name, 'complex'))
            #         # copy the xyz over from the parent dir - only one xyz in the dir, but no idea of the name- if _ in the name, the parent dir will be a number, or it might be the nsame of the complex? 
                    if 'equil.xyz' in listdir('.'):
                        system(f'cp equil.xyz {self.base_name}/complex/complex.xyz')
                system(f'mv {self.base_name}.inp {self.base_name}.job {self.base_name}/complex/')

    def make_run_dir(self):
        if not self.made_run_dir: # only do it once
            if not exists(self.base_name):
                mkdir(self.base_name) # make opt/spec/hessin parent dir
            self.made_run_dir = True


    def output_data(self, job_class, charge, mult, frags_in_subdir):
        self.create_inp()
        self.create_job()
        self.make_run_dir()
        self.place_files_in_dir()
        if frags_in_subdir:
            self.create_inputs_for_fragments(job_class, charge, mult)


    def create_inputs_for_fragments(self, job_class, charge, mult):
        """Very useful to generate files for each fragment automatically, for single point and frequency calculations, generating free energy changes. Called if ``frags_in_subdir`` is set to True, as each fragment is given a subdirectory in an overall subdirectory, creating the following directory structure (here for a 5-molecule system):
            .
            ├── frags
            │   ├── acetate0
            │   │   ├── acetate0.xyz
            │   │   └── spec.inp
            │   ├── acetate1
            │   │   ├── acetate1.xyz
            │   │   └── spec.inp
            │   ├── choline2
            │   │   ├── choline2.xyz
            │   │   └── spec.inp
            │   ├── choline3
            │   │   ├── choline3.xyz
            │   │   └── spec.inp
            │   └── water4
            │       ├── spec.inp
            │       └── water4.xyz
            ├── spec.inp
        """
        self.is_complex = False
        #look over self.mol.fragments, generate inputs- make a settings object with the desired features
        if not hasattr(self.mol, 'fragments'):
            self.mol.separate()
        #make subdir if not already there
        subdirectory = join(getcwd(), self.base_name, 'frags')
        if not exists(subdirectory):
            mkdir(subdirectory)

        parent_dir = getcwd()
        count = 0 #avoid  overwriting files by iterating with a number
        for frag, data in self.mol.fragments.items():
            if data['frag_type'] == 'frag':

                #make a directory inside the subdir for each fragment
                name = f"{data['name']}_{count}" # i.e. acetate0, acetate1, choline2, choline3, water4
                if not exists(join(subdirectory, name)):
                    mkdir(join(subdirectory, name)) # ./frags/water4/
                chdir(join(subdirectory, name))
                write_xyz(atoms = data['atoms'], filename = name + str('.xyz'))
            
                # re-use settings from complex
                if hasattr(self, 'merged'):
                    frag_settings = self.merged
                else:
                    frag_settings = self.defaults
                frag_settings.input.charge = data['charge']
                if data['multiplicity'] != 1:
                    frag_settings.input.multiplicity = data['multiplicity']
                job = job_class(using = name + str('.xyz'), settings=frag_settings, run_dir = True) 
                chdir(parent_dir)
                count += 1

        if hasattr(self.mol, 'ionic'):
            # only 1 ionic network        
            subdir_ionic = join(getcwd(), self.base_name, 'ionic')
            if not exists(subdir_ionic):
                mkdir(subdir_ionic)
            chdir(subdir_ionic)
            write_xyz(atoms = self.mol.ionic['atoms'], filename = 'ionic.xyz')
        
            # re-use settings from complex
            if hasattr(self, 'merged'):
                frag_settings = self.merged
            else:
                frag_settings = self.defaults
            frag_settings.input.charge = self.mol.ionic['charge']
            if self.mol.ionic['multiplicity'] != 1:
                frag_settings.input.multiplicity = self.mol.ionic['multiplicity']
            print('Creating input for the ionic network...')
            job = job_class(using = 'ionic.xyz', settings=frag_settings, fmo = True, run_dir = True) 
            chdir(parent_dir)

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
        elif self.sc == 'mon':
            jobfile = jobfile.replace('base_name', f'{self.base_name}') 
        self.write_file(jobfile, filetype='job')