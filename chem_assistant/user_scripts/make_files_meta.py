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
    parent = os.getcwd()
    for path, dirs, files in os.walk(base_dir):
        for file in files:
            if file == 'meta.py':
                os.chdir(path)
                if not any(file.endswith('.xyz') for file in os.listdir('.')):
                    raise TypeError(f'Meta file requires an xyz file in the same directory. Check {os.getcwd()}')
                else:
                    for file in os.listdir('.'):
                        if file.endswith('.xyz'):
                            print(os.path.join(os.getcwd(), file))
                os.system('python3 meta.py')
                os.chdir(parent)
        
            



