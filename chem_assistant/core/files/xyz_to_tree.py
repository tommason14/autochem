#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File: xyz_to_tree.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Script to take a directory of xyz files and create a directory structure with input and
job files 
"""
import os
from shutil import copyfile

def check_dir():
    """If in files directory, do nothing. If a subdir is called files, then move into that"""
    if os.getcwd().split('/')[-1] == 'files':
        pass
    elif os.path.isdir('files'):
        os.chdir(os.path.join(os.getcwd(), 'files'))

def get_xyz():
    return [file for file in os.listdir('.') if file.endswith('.xyz')]

def make_tree(files):
    #move up one dir, so you have a structure of file separate to the tree 
    # .
    # |
    # |__files
    # |__tree
    #       |
    #       |__ subdirs 
    os.chdir(os.path.dirname(os.getcwd())) #getting the parent dir of the files dir, then moving to it
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
        os.chdir(cwd)       

check_dir()
files = get_xyz()
make_tree(files)
#print(os.getcwd()) #now in 'calcs' dir
