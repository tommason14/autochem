#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File: psi.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Interface between Python and PSI4 input files 
"""

from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.settings import (Settings, read_template, dict_to_settings)
from ..core.job import Job
from ..core.periodic_table import PeriodicTable as PT
from ..core.sc import Supercomp

from os import (chdir, mkdir, getcwd)
from os.path import (exists, join, dirname)

__all__ = ['PsiJob']

class PsiJob(Job):
    # Note the job scripts require the supercomputer to be entered, such as:

    # >>> j = PsiJob(using = 'file.xyz')
    # >>> j.supercomp = 'raijin'
    """Class for creating PSI4 input files and job scripts. 
    
    The input files generated default to single point energy calculations using MP2/cc-pVTZ, with frozen core orbitals- for this reason, these single point calculations are very fast. This is easily changed by creating a |Settings| object and adding parameters, using the following syntax.

    >>> s = Settings()
    >>> s.input.globals.basis= 'cc-pVDZ'
    >>> s.input.molecule.extra_value = 'extra'
    >>> j = PsiJob(using = '../xyz_files/mesylate.xyz', settings = s)

    This yields the following result:

        memory 32 Gb

        molecule complex {
        -1 1
         C      -11.52615475      2.13587901     -3.92614475
         H      -12.17727298      1.37283268     -4.39314733
         H      -12.13111156      2.84527650     -3.33020803
         H      -10.95289836      2.67720525     -4.70258006
         S      -10.36648767      1.31567304     -2.82897636
         O       -9.54405868      2.38757303     -2.22205822
         O      -11.24567273      0.60890457     -1.83183396
         O       -9.60100212      0.36690604     -3.68579623
        units angstrom
        no_reorient
        symmetry c1
        extra_value extra
        }

        set globals {
            basis cc-pVDZ
            scf_type DF
            freeze_core True
            guess sad
            S_ORTHOGONALIZATION canonical
        }

        energy('mp2')

    Options are added in sections: 
        - self.input.molecule for any key value pair in the molecule section
        - self.input.globals for the 'set globals' section
        - any options not enclosed in braces appear before the last line
        - To change the run type:
            >>> self.input.run = {'optimize': 'scf'}
            # produces optimize('scf')
        - If extra run options are required:
            >>> self.input.run.additional = {'dertype': 'energy'} 
            # produces optimize('scf', dertype='energy')
        NOTE: An option for adding commands outside of the molecule and globals section needs to be added.
                                
    
    The names of files created default to the type of calculation: optimisation (opt), single point energy (spec) or hessian matrix calculation for thermochemical data and vibrational frequencies (freq). If a different name is desired, pass a string with the ``filename`` parameter, with no extension. The name will be used for both input and job files.
        >>> job = GamessJob(using = 'file.xyz', filename = 'benzene')
    This command produces two files, benzene.inp and benzene.job.
    
    If a system is comprised of multiple fragments, each fragment can have its own input file created in a subdirectory by passing in ``frags_in_subdir`` = True.
    """
    def __init__(self, using = None, frags_in_subdir = False, settings = None, filename = None):
        super().__init__(using)
        self.filename = filename
        self.defaults = read_template('psi.json') #settings object 
        if settings is not None:
            # can only have one run type, currently- need to delete the energy if running an
            # optimisation, for example
            if 'run' in settings.input.keys():
                del self.defaults.input.run
            self.merged = self.defaults.merge(settings) # merges inp, job data 
            self.input = self.merged.input
            self.job = self.merged.job
        else:
            self.input = self.defaults.input
        if '/' in using:
            self.title = using.split('/')[-1][:-4] #say using = ../xyz_files/file.xyz --> 
        else:
            self.title = using[:-4]
        
        self.create_inp()
        self.create_job()
        if frags_in_subdir:
            self.create_inputs_for_fragments()
        
    def make_header(self):
        """ Transform all contents of |Settings| objects into PSI4 input file headers, containing all the information pertinent to the calculation
         """
        mem  = f"memory {self.input.memory}\n\n"
        mol = "molecule complex {\n"
        charge = f"{self.input.molecule.charge} {self.input.molecule.multiplicity}\n"
        atoms = ""
        for atom in self.mol.coords:
            atoms += f" {atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}\n"
        units = f"units {self.input.molecule.units}\n"
        sym = f"symmetry {self.input.molecule.symmetry}\n"
        reorient = "no_reorient\n"
        end = "}\n"
                
        data = [mem, mol, charge, atoms, units, reorient, sym,  end]
        
        # add in user options
        for key, value in self.input.molecule.items():
            if key not in ("charge", "multiplicity", "units", "symmetry"):
                key = f"{key} {value}\n"
                data.insert(-1, key) #insert before last item
        self.header = data
       
    def add_globals(self):
        inp = self.header
        inp.append('\nset globals {\n')
        for key, value in self.input.globals.items():
            inp.append(f"    {key} {value}\n")
        inp.append('}\n')
        self.inp = inp
    
    def add_run(self):
        res = []
        # list of tuples- to ensure the 'normal' entry, the one defined input.run, appears first
        # in the list by adding a counter. In testing, {optimize: scf} with additional {dertype:
        # energy} produced dertype('energy', optimize='scf'), not optimize('scf', dertype='energy')
        # due to alphabetical ordering of [('dertype', 'energy'), ('optimize', 'scf')]
        # probably should have just made two lists of tuples, one for normal, one for additional
        for k,v in self.input.run.items():
            if k != 'additional':
                res.append((0, k, v))
            if k == 'additional': 
                counter = 1
                for k1, v1 in self.input.run[k].items():
                    res.append((counter, k1, v1)) 
                    counter += 1
                # if I ever need to add two different types of run in the same file,
                # comment out the two lines above, add in the 3 below, and add an additional
                # dict key: (AND CHANGE DOCSTRING)
                # would need to add resulting string to a list and then concatenate that list with
                # self.inp as well, otherwise you would have a combination in the same line
                ##############
                # s = Settings()
                # s.input.run = {'optimize': 'scf'}
                # s.input.run.additional = {'optimize' :{'dertype': 'energy',
                #                                        'entry': 'value'}}
                ##############
                # for data in self.input.run[k].values():
                #     for k1, v1 in data.items():
                #         res.append((k1, v1))
        res = sorted(res, key = lambda val: val[0]) #sort by the first item of tuple, the number
        string = f"{res[0][1]}('{res[0][2]}'"
        for val in res[1:]:
            string += f", {val[1]}='{val[2]}'"
        string += ')'
        self.inp.append(string)
        self.inp = "".join(self.inp)
 
    def write_file(self, data, filetype):
        """Writes the generated GAMESS input/jobs to a file. If no filename is passed when the class is instantiated, the name of the file defaults to the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (freq).""" 
        for key in self.input.run.keys(): #run, or additional
            if key != 'additional':
                nom = key
        if self.filename == None:
            options = {'optimize': 'opt', 'energy': 'spec', 'frequency': 'freq'}
            name = options.get(nom, 'file') #default name = file
        else:
            name = self.filename
        with open(f"{name}.{filetype}", "w") as f:
            f.write(data)

    def create_inp(self):
        self.make_header()
        self.add_globals()
        self.add_run()
        self.write_file(self.inp, filetype = 'inp')

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

    def create_job(self):
        """Returns the relevant job template as a list, then performs the necessary modifications. After, the job file is printed in the appropriate directory."""
        # RJN:
        # job memory = 
        # ncpus = nfrags
        # jobfs = 
        # MGS:

        # GAIA:

        self.get_sc()
        job = f"gamess_{self.sc}.job"
        dir_name = dirname(__file__) # interfaces (dir of gamess.py)
        templates = join(dir_name, '..', 'templates')
        job_file = join(templates, job)

        job = []
        with open(job_file) as f:
            for line in f:
                job.append(line)        

        # modify

        job = "".join(job)
        # write
        self.write_file(job, filetype="job")               

    def create_inputs_for_fragments(self):
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
       
        self.mol.separate() #creating frags 
        #look over self.mol.fragments, generate inputs- make a settings object with the desired features

        #make subdir if not already there
        subdirectory = join(getcwd(), 'frags')
        if not exists(subdirectory):
            mkdir(subdirectory)

        parent_dir = getcwd()
        count = 0 #avoid  overwriting files by iterating with a number
        for frag, data in self.mol.fragments.items():
            #make a directory inside the subdir for each fragment
            name = f"{data['name']}_{count}" # i.e. acetate0, acetate1, choline2, choline3, water4
            if not exists(join(subdirectory, name)):
                mkdir(join(subdirectory, name)) # ./frags/water4/
            chdir(join(subdirectory, name))
            Molecule.write_xyz(self, atoms = data['atoms'], filename = name + str('.xyz')) #using the method, but with no class
            
            #use the same settings, so if runtype is freq, generate freq inputs for all fragments too.
            if hasattr(self, 'merged'):
                frag_settings = self.merged
            else:
                frag_settings = self.defaults
            frag_settings.input.molecule.charge = data['charge']
            if data['multiplicity'] != 1:
                frag_settings.input.molecule.multiplicity = data['multiplicity']
            job = PsiJob(using = name + str('.xyz'), settings=frag_settings) 
            chdir(parent_dir)
            count += 1
            
        

