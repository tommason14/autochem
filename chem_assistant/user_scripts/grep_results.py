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
from ..core.results import (GamessResults, PsiResults)
import os


def get_logs(directory):
    logs = []
    for path, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.log') or file.endswith('.out'):
                logs.append(os.path.join(path, file))
    return logs


def make_instance(log):
    """Return an instance of the desired class- |GamessResults|, |PsiResults|"""
    log_type = get_type(log)
    logs = {'gamess': GamessResults(log),
            'psi4': PsiResults(log)}
    return logs[log_type]


def fetch_data(log):
    r = make_instance(log)
    if r.completed():
        if r.get_runtype() == 'optimize':
            r.get_equil_coords()
        en = r.get_energy() #fmo3 > fmo2 > non-fmo #fix gamess energy, running v slow
        return log, en
    else:
        r.get_error()


def parse_results(dir):
    cwd = os.getcwd()
    for log in get_logs(dir):
        f = fetch_data(log)
        print(f)

            
