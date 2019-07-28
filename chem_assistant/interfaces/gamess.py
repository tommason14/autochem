#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File: gamess.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Interface between Python and creating GAMESS input files 
"""

from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.settings import (Settings, read_template, dict_to_settings)
from ..core.job import Job
from ..core.periodic_table import PeriodicTable as PT
from ..core.sc import Supercomp
from ..core.utils import (sort_elements, write_xyz)

from os import (chdir, mkdir, getcwd, system, walk, listdir)
from os.path import (exists, join, dirname)

__all__ = ['GamessJob']

class GamessJob(Job):
    # Note the job scripts require the supercomputer to be entered, such as:

    # >>> j = GamessJob(using = 'file.xyz')
    # >>> j.supercomp = 'raijin'
    """Class for creating GAMESS input files and job scripts. 
    
    The input files generated default to geometry optimisations at the SRS-MP2/cc-pVDZ level of theory. This is easily changed by creating a |Settings| object and adding parameters, using the following syntax.
    >>> s = Settings()
    >>> s.input.contrl.runtyp = 'energy'
    >>> s.input.basis.gbasis = 'CCT'
    >>> s.input.mp2.scsopo = 1.64
    >>> j = GamessJob(using = '../xyz_files/ch_ac.xyz', settings = s)
    This yields the following result:
         $SYSTEM MEMDDI=0 MWORDS=500 $END
         $CONTRL ICHARG=0 ISPHER=1 MAXIT=200 RUNTYP=ENERGY SCFTYP=RHF $END
         $STATPT NSTEP=500 $END
         $SCF DIIS=.TRUE. DIRSCF=.TRUE. FDIFF=.FALSE. $END
         $BASIS GBASIS=CCT $END
         $MP2 CODE=IMS SCSOPO=1.64 SCSPAR=0.0 SCSPT=SCS $END
         $DATA
        ch_ac
        C1
         H 1.0
         C 6.0
         N 7.0
         O 8.0
         $END
         C     6.0   -6.27719   -2.98190   -7.44828
         H     1.0   -7.26422   -2.51909   -7.63925
         H     1.0   -5.81544   -3.32313   -8.39516 
    If FMO (Fragment Molecular Orbital) calculations are desired, pass the keyword argument ``fmo``, set to *True*, along with a settings object with a parameter of ``nfrags``:
        >>> s = Settings()
        >>> s.nfrags = 4
        >>> job = GamessJob(using = 'file.xyz', fmo = True, settings = s)
    
    The class creates different subdirectories for every molecule in the system.
    Using the class with `frags_in_subdir` set to true produces:
        - a `complex` subdirectory- one per xyz
        - an `ionic` subdirectory for every complex with the neutral species and/or single atom ions removed
            - for water inclusion, N2 inclusion, alkali metal inclusion
        - a `frags` subdirectory for every fragment of the complex
    
    The names of files created default to the type of calculation: optimisation (opt), single point
energy (spec) or hessian matrix calculation for thermochemical data and vibrational frequencies
(hess). If a different name is desired, pass a string with the ``filename`` parameter, with no extension. The name will be used for both input and job files.
        >>> job = GamessJob(using = 'file.xyz', fmo = True, filename = 'benzene')
    
    This command produces two files, benzene.inp and benzene.job.
    
    Files are placed in a subdirectory of their own name. So when creating optimisation files, files are placed in opt:
        .
        └── opt
            ├── opt.inp
            └── opt.job
    """
    def __init__(self, using = None, fmo = False, frags_in_subdir = False, settings = None, filename
= None, is_complex = False, run_dir = None):
        super().__init__(using)
        self.fmo = fmo # Boolean
        self.filename = filename
        self.defaults = read_template('gamess.json') #settings object 
        if settings is not None:
            self.merged = self.defaults.merge(settings) # merges inp, job data 
            self.input = self.merged.input
            self.job = self.merged.job
        else:
            self.input = self.defaults.input
        if '/' in using:
            self.title = using.split('/')[-1][:-4] #say using = ../xyz_files/file.xyz --> 
        else:
            self.title = using[:-4]
        
        self.is_complex = is_complex # creates a `complex` dir

        if run_dir is not None:
            self.made_run_dir = True
        else:
            self.made_run_dir = False
         
        self.create_inp()
        self.create_job()
        if frags_in_subdir:
            self.create_inputs_for_fragments(complex_is_fmo = self.fmo)
        
    def determine_fragments(self):
        if self.fmo:
            self.mol.separate()
            fmo_data = self.fmo_formatting()
            self.input.fmo = fmo_data
            self.input.fmoprp.maxit = 200
            # num_ions = len([frag['name'] for frag in self.mol.fragments.values() if frag['name'] not in Molecule.Neutrals])
            # self.input.gddi.ngroup = num_ions # for hy2ip use 4 groups not 5
            self.input.gddi.ngroup = len(self.mol.fragments)

    def order_header(self):
        if self.fmo:
            desired = ['SYSTEM', 'CONTRL', 'GDDI', 'STATPT', 'SCF', 'BASIS', 'FMO', 'FMOPRP'] # mp2/dft after
        else:
            desired = ['SYSTEM', 'CONTRL', 'STATPT', 'SCF', 'BASIS']


        self.header = []
        for i in desired:
            for line in self.unordered_header:
                item = line.split()[0][1:]
                if item == i:
                    self.header.append(line)
        # add in additional commands after, given as sett.input.blah = 'blah'
        else:
            for line in self.unordered_header:
                item = line.split()[0][1:]
                if item not in desired:
                    self.header.append(line)
        
        # no need for stationary point steps if not an optimisation
        if self.input.contrl.runtyp != 'optimize':
            for index, line in enumerate(self.header):
                if line.split()[0]  == '$STATPT':
                    del self.header[index]

        self.header = "".join(self.header)

    def make_automatic_changes(self):
        """
        Common scenarios are implemented to remove the commands needed to be called by the user.
        For example, this sets the opposite spin parameter for SRS-MP2 for commonly used correlation
        consistent basis sets without the need to define it in the user code.
        """
        opp_spin_params = {'cct': 1.64,
                           'ccq': 1.689,
                           'accd': 1.372,
                           'acct': 1.443,
                           'accq': 1.591}
        if 'mp2' in self.input:
            for basis, opp in opp_spin_params.items():
                if self.input.basis.gbasis.lower() == basis:
                    self.input.mp2.scsopo = opp
                    break
        if 'fmo' and 'pcm' in self.input:
            self.input.pcm.ifmo = -1

    def parse_settings(self):
        """Transforms all contents of |Settings| objects into GAMESS input file headers, containing all the information pertinent to the calculation"""
        
        def format_line_if_too_long(line):
            if len(line) > 55:
                line = line.split()
                # find suitable character length to split on
                line_length = 0
                chars = []
                for val in line:
                    line_length += len(val)
                    chars.append(line_length)
                for index, length in enumerate(chars):
                    if length > 55:
                        line.insert(index - 1, '\n ')
                        break # only insert one newline, and indent
                # space before, newline at end
                line.insert(0, '')
                line = ' '.join(line) + '\n'
            return line

        def parse(key, value):
            ret = ''
            if isinstance(value, Settings):
                ret += ' ${}'.format(key.upper())
                for el in value:
                    ret += ' {}={}'.format(el.upper(), str(value[el]).upper())
                ret += ' $END\n'
            else:
                ret += ' ${} {}\n $END\n'.format(key.upper(), value.upper())
            if key is not 'fmo':
                ret = format_line_if_too_long(ret)
            return ret
    

        inp = [parse(item, self.input[item])
            for item in self.input]
        return inp

    #format for fmo run of complete system 
    def fmo_meta(self):
        """Creates strings for the INDAT and ICHARG blocks of GAMESS FMO calculations, bound to the
        molecule instance as self.indat and self.charg"""
        info = {}
        #group together indat and charge for each fragment, and order according to the atom indices of indat 
        for frag, data in self.mol.fragments.items():
            if frag is not 'ionic':
                if len(data['atoms']) == 1:
                    info[frag] = {"indat": f"0,{data['atoms'][0].index},-{data['atoms'][0].index},",
                    "charg" : str(data['charge'])}
                else:
                    info[frag] = {"indat": f"0,{data['atoms'][0].index},-{data['atoms'][-1].index},",
                    "charg" : str(data['charge'])}
        # items need sorting
        # 0,1,7, ### sort on 2nd item ###
        # 0,8,28,
        # 0,29,35,
        # and not 
        # 0,1,7,
        # 0,29,35,
        # 0.8,28 
        
        sorted_info = sorted(info.items(), key = lambda val: int(val[1]['indat'].split(',')[1])) 
        # could also just sort on mol (or frag of info[frag]), from the assignments in self.split(), but these might not
        # always be in a numerical order- by using the index from self.coords, it is always ensured that the
        # correct order is shown, as these coords are also used in the input file        
        self.indat = []
        self.charg = []
        for val in sorted_info:
            self.indat.append(val[1]['indat'])
            self.charg.append(val[1]['charg'])
        self.indat.append('0')

    def fmo_formatting(self):
        self.fmo_meta() # gives self.mol.indat, self.mol.charg
        if self.input.contrl.runtyp.lower() in ('optimize', 'hessian', 'fmohess'):
            nbody = 2
            rcorsd = 100
        else:
            #FMO3 for specs- change rcorsd to 50
            nbody= 3
            rcorsd = 50
        
        string = f"\n     NFRAG={len(self.mol.fragments)} NBODY={nbody}\n"
        if 'mp2' in self.input:
            string += "     MPLEVL(1)=2\n"
        string += f"     INDAT(1)={self.indat[0]}\n"
        for d in self.indat[1:]:
            string += f"{' '*14}{d}\n"
        string += f"     ICHARG(1)={','.join(self.charg)}\n"
        string += f"     RESPAP=0 RESPPC=-1 RESDIM=100 RCORSD={rcorsd}"
        return string

    def make_inp(self):
        inp = self.header
        inp += ' $DATA\n'
        inp += f'{self.title}\n'
        inp += 'C1\n'
        if self.fmo:
            for el in self.mol.complex['elements']: #list of tuples [('H', 1.0), ('O', 8.0)]
                inp += f" {el[0]} {el[1]}\n"
            inp += " $END\n"
            inp += " $FMOXYZ\n"
        for atom in self.mol.coords:
            inp += f" {atom.symbol:5s} {PT.get_atnum(atom):>3}.0{atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}\n"
        inp += ' $END'
        return inp

    def file_basename(self):
        """If no filename is passed when the class is instantiated, the name of the file defaults to
        the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (hess). This method creates an attribute ``base_name``, used in creating the input and job files."""

        if self.filename is not None:
            self.base_name = self.filename
        else:
            options = {'optimize': 'opt', 'energy': 'spec', 'hessian': 'hess', 'fmohess': 'hess'}
            self.base_name = options.get(self.input.contrl.runtyp, 'file') #default name = file

    def write_file(self, data, filetype):
        """Writes the generated GAMESS input/jobs to a file. If no filename is passed when the class is instantiated, the name of the file defaults to the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (freq). 
        NOTE: Must pass data as a string, not a list!""" 
        with open(f"{self.base_name}.{filetype}", "w") as f:
            f.write(data)

    def create_inp(self):
        self.input = self.input.remove_none_values()
        self.determine_fragments() #add fmo info to input settings, if self.fmo is True
        self.make_automatic_changes()
        self.unordered_header = self.parse_settings() 
        self.order_header() #create self.header variable
        inp = self.make_inp()
        self.file_basename()
        self.write_file(inp, filetype = 'inp')
   
    def get_job_template(self):
        job_file = self.find_job()
        with open(job_file) as f:
            job = f.read()       
            return job
    
    def change_mgs_job(self, job):
        if hasattr(self.mol, 'fragments') and len(self.mol.fragments) != 0:
            num_frags = len(self.mol.fragments)
            jobfile = job.replace('nodes=1', f'nodes={num_frags}')
            jobfile = jobfile.replace('24 24', f'{24 * num_frags} 24') 
            return jobfile
        return job
            
    def change_rjn_job(self, job):
        if hasattr(self.mol, 'fragments') and len(self.mol.fragments) != 0:
            num_frags = len(self.mol.fragments)
            jobfile = job.replace('ncpus=32', f'ncpus={16 * num_frags}')
            jobfile = jobfile.replace('mem=125gb', f'mem={4 * 16 * num_frags}gb') # 4gb cpus
            jobfile = jobfile.replace('jobfs=150gb', f'jobfs={4 * 16 * num_frags + 20}gb')
            return jobfile
        return job
    
    def create_job(self):
        """Returns the relevant job template as a list, then performs the necessary modifications. After, the job file is printed in the appropriate directory."""
        jobfile = self.get_job_template()
        if self.sc == 'mgs':
            jobfile = self.change_mgs_job(jobfile)
            jobfile = jobfile.replace('name', f'{self.base_name}') 
        elif self.sc == 'rjn':
            jobfile = self.change_rjn_job(jobfile)
            jobfile = jobfile.replace('name', f'{self.base_name}') 
        elif self.sc == 'mon':
            jobfile = jobfile.replace('base_name', f'{self.base_name}') 
        elif self.sc == 'mas':
            jobfile = jobfile.replace('base_name', f'{self.base_name}') 
        elif self.sc == 'stm':
            jobfile = jobfile.replace('name', f'{self.base_name}') 
        self.write_file(jobfile, filetype='job')

    def make_run_dir(self):
        if not self.made_run_dir: # only do it once
            if not exists(self.base_name):
                mkdir(self.base_name) # make opt/spec/hessin parent dir
            self.made_run_dir = True

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

    def create_inputs_for_fragments(self, complex_is_fmo = False):
        """Very useful to generate files for each fragment automatically, for single point and frequency calculations, generating free energy changes. Called if ``frags_in_subdir`` is set to True, as each fragment is given a subdirectory in an overall subdirectory, creating the following directory structure (here for a 5-molecule system):
            .
            ├── frags
            │   ├── acetate0
            │   │   ├── acetate0.xyz
            │   │   └── spec.inp
            │   ├── acetate1
            │   │   ├── acetate1.xyz
            │   │   └── spec.inp
            │   ├── choline2
            │   │   ├── choline2.xyz
            │   │   └── spec.inp
            │   ├── choline3
            │   │   ├── choline3.xyz
            │   │   └── spec.inp
            │   └── water4
            │       ├── spec.inp
            │       └── water4.xyz
            ├── spec.inp
        """
        self.is_complex = False
        #look over self.mol.fragments, generate inputs- make a settings object with the desired features
        if not hasattr(self.mol, 'fragments'):
            self.mol.separate()
        #make subdir if not already there
        subdirectory = join(getcwd(), 'frags')
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
                frag_settings.input.contrl.icharg = data['charge']
                if data['multiplicity'] != 1:
                    frag_settings.input.contrl.mult = data['multiplicity']
                job = GamessJob(using = name + str('.xyz'), settings=frag_settings, run_dir = True) 
                chdir(parent_dir)
                count += 1

        if hasattr(self.mol, 'ionic'):
            # only 1 ionic network        
            subdir_ionic = join(getcwd(), 'ionic')
            if not exists(subdir_ionic):
                mkdir(subdir_ionic)
            chdir(subdir_ionic)
            write_xyz(atoms = self.mol.ionic['atoms'], filename = 'ionic.xyz')
        
            # re-use settings from complex
            if hasattr(self, 'merged'):
                frag_settings = self.merged
            else:
                frag_settings = self.defaults
            frag_settings.input.contrl.icharg = self.mol.ionic['charge']
            if self.mol.ionic['multiplicity'] != 1:
                frag_settings.input.contrl.mult = self.mol.ionic['multiplicity']
            ## FMO only if more than 2 fragments
            
            print(Molecule(atoms=self.mol.ionic['atoms']).separate())
            job = GamessJob(using = 'ionic.xyz', settings=frag_settings, fmo = complex_is_fmo, run_dir = True) 
            chdir(parent_dir)
