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
import os
from shutil import copyfile

# from gamess import GamessJob
# from psi import PsiJob

def check_dir():
    """If in files directory, do nothing. If a subdir is called files, then move into that"""
    if os.getcwd().split('/')[-1] == 'files':
        pass
    elif os.path.isdir('files'):
        os.chdir(os.path.join(os.getcwd(), 'files'))

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


def make_tree(files):
    #move up one dir, so you have a structure of file separate to the tree 
    # .
    # |
    # |__files
    # |__tree
    #       |
    #       |__ subdirs 
    # os.chdir(os.path.dirname(os.getcwd())) #getting the parent dir of the files dir, then moving to it
    calc_dir = os.path.join(os.getcwd(), 'calcs')
    if not os.path.exists(calc_dir):
        os.mkdir(calc_dir)
    os.chdir(calc_dir)
    cwd = os.getcwd()
    for file in files:
        filename = file[:-4] # rm .xyz 
        new_dirs = []
        # store each portion in a list, then iterate over the list making dirs if needed
        for part in filename.split('_'):
            new_dirs.append(part)
        for idx, d in enumerate(new_dirs):
            new_dir = os.path.join(os.getcwd(), d)
            if os.path.isdir(d):
                os.chdir(new_dir)
            elif not os.path.isdir(d):
                os.mkdir(new_dir)
                os.chdir(new_dir)
            if idx + 1  == len(new_dirs): #when at maximum depth
                xyz = os.path.join(os.path.dirname(cwd), file)
                dest = os.path.join(os.getcwd(), file)
                copyfile(xyz, dest)
                # make_files(choice, file, settings) #calls GamessJob etc....
        os.chdir(cwd)
       
def xyz_to_tree(settings = None):
    # job_choice = get_choices()
    check_dir()
    files = get_xyz()
    make_tree(files)

# xyz_to_tree()
