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
