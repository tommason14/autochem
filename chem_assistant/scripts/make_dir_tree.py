#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File: xyz_to_tree.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Script to take a directory of xyz files and create a directory structure with input and
job files. To use, need to import module, and add settings:
import xyz_to_tree
s = Settings()
s.....

xyz_to_tree(settings = s) 
"""
from ..interfaces.gamess import GamessJob
# Add user inputs- packages? fmo? 

import os
from shutil import copyfile

__all__ = ['xyz_to_tree']

def check_dir():
    """If in files directory, do nothing. If a subdir is called files, then move into that"""
    if os.getcwd().split('/')[-1] == 'files':
        return os.getcwd()
    elif os.path.isdir('files'):
        os.chdir(os.path.join(os.getcwd(), 'files'))
        return os.getcwd()

def get_xyz():
    return [file for file in os.listdir('.') if file.endswith('.xyz')]

def ask_package():
    print('Which software would you like?')
    print("""\
1. GAMESS
2. PSI4
3. LAMMPS""")
    done = False
    while not done:
        choice = int(input('Choice [1,2,3]: '))
        if choice in (1,2,3):
            done = True
        else:
            print('Please choose 1-3')
    return choice    

def fmo_for_gamess():
    ret = input('Run FMO calculations? [Y/N] ')
    if ret.lower() not in ('y', 'n'):
        print('Please choose Y or N')
        fmo_for_gamess()
    return ret.lower()

def get_choices():
    c = ask_package()
    if c == 1:
        fmo = fmo_for_gamess()
        if fmo == 'y':
            return 'gamess_fmo'
        else:
            return 'gamess'
    elif c == 2:
        return 'psi'
    elif c == 3:
        return 'lammps'

def make_files(choice, xyz, s):
    choices = {
        "gamess": GamessJob(using = xyz, frags_in_subdir = True, settings = s),
        "gamess_fmo": GamessJob(using = xyz, fmo = True, frags_in_subdir = True, settings = s),
        "psi": PsiJob(using = xyz, frags_in_subdir = True, settings = s),
        "lammps": None
    }    
    return choices[choice] 


def make_parent_dir():
    calc_dir = os.path.join(os.path.dirname(os.getcwd()), 'calcs')
    if not os.path.exists(calc_dir):
        os.mkdir(calc_dir)
    os.chdir(calc_dir)
    return os.getcwd()

def make_dir_list(file):
    filename = file[:-4] # rm .xyz 
    new_dirs = []
        # store each portion in a list, then iterate over the list making dirs if needed
    for part in filename.split('_'):
        new_dirs.append(part)
    return new_dirs

def change_to_subdir(subdirectory):
    new_dir = os.path.join(os.getcwd(), subdirectory)
    if os.path.isdir(new_dir):
        os.chdir(new_dir)
    elif not os.path.isdir(new_dir):
        os.mkdir(new_dir)
        os.chdir(new_dir)

def copy_xyz(xyz_dir, file):
    xyz = os.path.join(xyz_dir, file)
    dest = os.path.join(os.getcwd(), file)
    copyfile(xyz, dest)

def make_tree_and_copy(xyz_dir, files):
    """Makes a sibling directory to 'files' (named 'calcs'), creates subdirectories with names based on the xyz files in 'files', and then copies the xyz files from 'files' to the deepest sub directory of the path created.

    The end result is this:
    .
    ├── calcs
    │   ├── c1mim
    │   │   └── nh3
    │   │       ├── c1mim_nh3.xyz
    │   │       └── frags
    │   │           ├── c1mim_0
    │   │           │   └── c1mim_0.xyz
    │   │           ├── nh3_1
    │   │           │   └── nh3_1.xyz
    │   │           ├── nh3_2
    │   │           │   └── nh3_2.xyz
    │   │           └── nh3_3
    │   │               └── nh3_3.xyz
    │   ├── ch
    │   │   └── ac
    │   │       ├── nowater
    │   │       │   ├── ch_ac_nowater.xyz
    │   │       │   └── frags
    │   │       │       ├── acetate_0
    │   │       │       │   └── acetate_0.xyz
    │   │       │       ├── acetate_1
    │   │       │       │   └── acetate_1.xyz
    │   │       │       ├── choline_2
    │   │       │       │   └── choline_2.xyz
    │   │       │       └── choline_3
    │   │       │           └── choline_3.xyz
    │   │       └── water
    │   │           ├── ch_ac_water.xyz
    │   │           └── frags
    │   │               ├── acetate_0
    │   │               │   └── acetate_0.xyz
    │   │               ├── acetate_1
    │   │               │   └── acetate_1.xyz
    │   │               ├── choline_2
    │   │               │   └── choline_2.xyz
    │   │               ├── choline_3
    │   │               │   └── choline_3.xyz
    │   │               └── water_4
    │   │                   └── water_4.xyz
    │   └── water
    │       ├── frags
    │       │   ├── water_0
    │       │   │   └── water_0.xyz
    │       │   └── water_1
    │       │       └── water_1.xyz
    │       └── water.xyz
    └── files
        ├── c1mim_nh3.xyz
        ├── ch_ac_nowater.xyz
        ├── ch_ac_water.xyz
        └── water.xyz"""

    calc_dir = make_parent_dir()
    cwd = os.getcwd()
    for file in files:
        new_dirs = make_dir_list(file)
        for idx, d in enumerate(new_dirs):
            change_to_subdir(d)
            if idx + 1  == len(new_dirs): #when at maximum depth
                copy_xyz(xyz_dir, file)
                # make_files(choice, file, settings) #calls GamessJob etc....
        os.chdir(cwd)
    return calc_dir
    
def make_job_files(base_dir):
    parent = os.getcwd()
    for path, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.xyz'):
                os.chdir(path)
                print(file)
                GamessJob(using=file, fmo=True, frags_in_subdir = True)
                os.chdir(parent)
    
def xyz_to_tree(settings = None):
    # job_choice = get_choices()
    xyz_directory = check_dir()
    files = get_xyz()
    calc_dir = make_tree_and_copy(xyz_directory, files)
    make_job_files(calc_dir)
    os.chdir(xyz_directory)