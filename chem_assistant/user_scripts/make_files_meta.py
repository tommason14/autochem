#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File: make_files_meta.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Make inputs using a meta.py file in a directory with xyz files 
"""

__all__ = ['make_files_from_meta']

import os
from subprocess import Popen

def make_files_from_meta(base_dir):
    """This function looks for ``meta.py`` in any subdirectory below the directory passed as an argument. ``meta.py`` is then run as a python file, with the idea of including data for a job. As long as there is an xyz file in the same directory as ``meta.py``,  the desired result is passed. 
    
    Usage:
        Files in directory:
            - ch_ac.xyz
            - meta.py
        
        meta.py:

            from chem_assistant import GamessJob
            GamessJob(using='ch_ac.xyz', frags_in_subdir = True)

    Resulting directory structure:
    
        dir
        ├── ch_ac.xyz
        ├── frags
        │   ├── acetate_0
        │   │   ├── acetate_0.xyz
        │   │   └── opt
        │   │       ├── opt.inp
        │   │       └── opt.job
        │   ├── acetate_1
        │   │   ├── acetate_1.xyz
        │   │   └── opt
        │   │       ├── opt.inp
        │   │       └── opt.job
        │   ├── choline_2
        │   │   ├── choline_2.xyz
        │   │   └── opt
        │   │       ├── opt.inp
        │   │       └── opt.job
        │   ├── choline_3
        │   │   ├── choline_3.xyz
        │   │   └── opt
        │   │       ├── opt.inp
        │   │       └── opt.job
        │   └── water_4
        │       ├── opt
        │       │   ├── opt.inp
        │       │   └── opt.job
        │       └── water_4.xyz
        ├── meta.py
        └── opt
            ├── opt.inp
            └── opt.job
    """
    parent = os.getcwd()
    for path, dirs, files in os.walk(base_dir):
        for file in files:
            if file == 'meta.py':
                os.chdir(path)
                if not any(file.endswith('.xyz') for file in os.listdir('.')):
                    raise TypeError(f'Meta file requires an xyz file in the same directory. Check {os.getcwd()}')
                else:
                    print(os.getcwd())
                os.system('python3 meta.py')
                os.chdir(parent)
        
            



