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
from ..interfaces.gamess_results import GamessResults
from ..interfaces.psi_results import PsiResults
import os


def get_logs(directory):
    logs = []
    for path, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.log') or file.endswith('.out') and file != 'freq.out': # freq.out used for thermo calculations with the fortran code
                logs.append(os.path.join(path, file))
    return logs


def make_instance(log):
    """Return an instance of the desired class- |GamessResults|, |PsiResults|"""
    log_type = get_type(log)
    logs = {'gamess': GamessResults(log),
            'psi4': PsiResults(log)}
    return logs.get(log_type, None)


def srs_output(r):
    """Returns parameters of SRS-MP2 calculations; HF energy, Opposite and Same spin parameters, as
    well as the overall SRS-MP2 energy, as a tuple""" 
    hf = r.get_hf()
    opp = r.get_e_os()
    same = r.get_e_ss()
    srs = r.get_srs()
    return hf, opp, same, srs

def fetch_energies(r):
    """Returns energies of both Gamess and Psi4 calculations. 
    If possible, a tuple is returned of the form:
        log file, basis set, HF energy, MP2 opposite spin energy, MP2 same spin energy, SRS-MP2
energy
    If that is not possible:
        log file, basis set, total energy"""
    if r.is_optimisation():
        r.get_equil_coords()
    if type(r) == GamessResults:
        if r.get_fmo_level() != 0:
            energy = srs_output(r)
        else:
           energy = r.get_non_fmo() # Total Energy = .... 
    if type(r) == PsiResults:
        energy = srs_output(r)
    return energy

def process_hessian(r):
   """Takes a hessian log file and produces thermochemical data. In future, extend to visualising
normal modes and plotting IR spectra using node intensities."""
    # generate freq.out, fort.10, then clean at end by removing from dir
    

def parse_results(dir):
    cwd = os.getcwd()
    for log in get_logs(dir):
        r = make_instance(log)
        if r.completed():
            basis = r.get_basis()
            if not r.is_hessian():
                print(f'{r.log}: Finding energies')
                f = fetch_energies(r)
                # adding to csv/df/database
                # print((log, basis) + f)
            else:
                print(f'{r.log}: Hessian calc')
                process_hessian(r)
        else:
            # incomplete cases
            r.get_error()

            
