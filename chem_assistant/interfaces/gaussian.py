#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File: gaussian.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Interface between Python and GAUSSIAN input files 
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

__all__ = ['GaussJob']

class GaussJob(Job):
    # Note the job scripts require the supercomputer to be entered, such as:

    # >>> j = GaussJob(using = 'file.xyz')
    # >>> j.supercomp = 'raijin'
    """Class for creating GAMESS input files and job scripts. 
    
    The class creates different subdirectories for every molecule in the system.
    Using the class with `frags_in_subdir` set to true produces:
        - a `complex` subdirectory- one per xyz
        - an `ionic` subdirectory for every complex with the neutral species and/or single atom ions removed
            - for water inclusion, N2 inclusion, alkali metal inclusion
        - a `frags` subdirectory for every fragment of the complex
    
    The names of files created default to the type of calculation: 
    optimisation (opt), single point energy (spec) or hessian matrix 
    calculation for thermochemical data and vibrational frequencies (hess). 
    If a different name is desired, pass a string with the ``filename`` 
    parameter, with no extension. The name will be used for both input and job
    files.
    
        >>> job = GamessJob(using = 'file.xyz', fmo = True, filename = 'benzene')
    
    This command produces two files, benzene.inp and benzene.job.
    
    """
    # sett.proc
    # sett.mem
    # sett....
    # sett.input.opt.scf = 'tight'
    # sett.input.freq
    # sett.input.td.nstates = 10
    # sett.input.td.root = 7
    # sett.input.grid = 'ultrafine'
    # sett.input.method = 'wB97xD'
    # sett.input.basis = 'aug-ccpvdz' # -> wB97xD/aug-ccpVDZ 

    def __init__(self, using = None, fmo = False, frags_in_subdir = False, settings = None, filename = None, is_complex = False, run_dir = None):
        super().__init__(using)
        self.filename = filename
        self.defaults = read_template('gauss.json') #settings object 
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
        self.xyz = using
        
        self.create_complex_dir_if_required(is_complex, frags_in_subdir)

        if run_dir is not None:
            self.made_run_dir = True
        else:
            self.made_run_dir = False
         
        self.create_inp()
        self.create_job()
        self.place_files_in_dir()

        if frags_in_subdir:
            self.create_inputs_for_fragments(complex_is_fmo = self.fmo)
        
    def create_complex_dir_if_required(self, is_complex, make_frags):
        self.is_complex = is_complex 
        if make_frags and not is_complex:
            self.is_complex = True