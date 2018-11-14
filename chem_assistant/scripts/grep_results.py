#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['parse_results']

"""
File: grep_results.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Searches all sub dirs for results 
"""
"""
- Create opt dir or spec or hess-  when creating the jobs (make another dir)
- If geom_opt has completed, extract equilibrium coords: save in opt folder as equil.xyz & if no
  spec directory available, create it. If no equil.xyz in spec dir, then copy equil.xyz over to spec dir, and make files
accordingly.

Get results.
...
...
...
...
Create single points for completed geom_opts? [y/n]
"""

from ..core.utils import (read_file, get_type)
import os


def get_logs(directory):
    logs = []
    parent = os.getcwd()
    for path, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.log') or file.endswith('.out'):
                logs.append(os.path.join(path, file)) # path is fine as the function uses a user-defined directory
    return logs


def parse_logs(directory):
    cwd = os.getcwd()
    for log in get_logs(directory):
        print(os.path.split(log))

def parse_results(dir):
    parse_logs(dir)

            # os.chdir(path)
            # log_type = get_type(file)
            # if log_type == 'gamess':
            #     r = GamessResults(file)
            # elif log_type == 'psi4':
            #
            #     r = PsiResults(file)
            #
            # if r.get_runtype() == 'optimize':
            #     r.get_equil_coords()
             
