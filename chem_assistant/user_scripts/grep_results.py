#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['results_table', 'parse_results', 'thermochemistry', 'get_results_class',
'search_for_coords']

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
from ..core.thermo import thermo_data
from ..core.utils import (read_file, get_type, get_files)
from ..interfaces.gamess_results import GamessResults
from ..interfaces.psi_results import PsiResults
import os
import pandas as pd


def get_results_class(log):
    """Return an instance of the desired class- |GamessResults|, |PsiResults|"""
    log_type = get_type(log)
    logs = {'gamess': GamessResults(log),
            'psi4': PsiResults(log)}
    return log_type, logs.get(log_type, None)

def search_for_coords(dir):
    for log in get_files(dir, ('.log', '.out')):
        _, r = get_results_class(log)
        if r.completed():
            if r.is_optimisation():
                print(f'{r.log}: Finding equilibrium coordinates...')
                r.get_equil_coords()

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
    orig = os.getcwd()
    if type(r) == GamessResults:
        if r.get_fmo_level() != 0:
            energy = srs_output(r)
        else:
            energy = (r.get_non_fmo(),) # Total Energy = .... 
    if type(r) == PsiResults:
        energy = srs_output(r)
    return energy


def parse_results(dir):
    output = []
    cwd = os.getcwd()
    for log in get_files(dir, ('.out', '.log')):
        log = log[2:] # no ./file, just file
        filetype, r = get_results_class(log)
        try:
            if r.completed():
                basis = r.get_basis()
                if not r.is_hessian():
                    print(f'{r.log}: Finding energies')
                    f = fetch_energies(r)
                    # adding to csv/df/database
                    output.append((log, filetype, basis) + f)
                else:
                    print(f'{r.log}: Hessian calc')
            else:
                # incomplete cases
                r.get_error()
        except AttributeError:
            continue
    return output

def make_table(results):
    '''Create a dataframe of energies and other relevant data extracted from each log file'''
    import pandas as pd
    files = []
    types = []
    sets = []
    hf_list = []
    opp_list = []
    same_list = []
    energies = []
    for res in results:
        filename, filetype, basis, *energy = res
        if len(energy) > 1:
            hf, opp, same, en = energy
        else:
            hf = 'NA'
            opp = 'NA'
            same = 'NA'
            en = energy[0]
        #path, file = os.path.split(filename)
        #if path is not '.':
        #    file = 'subdir/' + file
        files.append(filename)
        types.append(filetype)
        sets.append(basis)
        hf_list.append(hf)
        opp_list.append(opp)
        same_list.append(same)
        energies.append(en)
    
    df = pd.DataFrame({'File': files, 'Type': types, 'Basis': sets,
                       'HF': hf_list, 'Opp Spin': opp_list, 'Same Spin': same_list,
                       'SRS (if poss)': energies})
    return df

def results_table(dir):
    results = parse_results(dir)
    table = make_table(results)
    return table

def thermochemistry(dir):
    """Returns thermochemical data for all the relevant hessian log files in the given directory and
    subdirectories. Saves to excel file, thermo.xlsx
    Usage:
        >>> thermochemistry('.')
    """
    collected = \
    {
        'File': [],
        'ZPVE': [],
        'TC': [],
        'S elec': [],
        'S tran': [],
        'S rot': [],
        'S vib': [],
        'S tot': [],
        'TC - TS': []
    }

    for log in get_files(dir, ('.log', '.out')):
        _, r = get_results_class(log)
        if r.completed():
            if r.is_hessian():
                print(f'Thermo data for {log}')
                res = thermo_data(r.log)
                res['File'] = os.getcwd() +  log[1:]
                for k, v in res.items():
                    collected[k].append(v)

    # add units to dict keys   

    kj = ('ZPVE', 'TC', 'TC - TS')
    jmol = ('S elec', 'S tran', 'S rot', 'S vib', 'S tot')
    kv = list(collected.items())
    collected.clear() 
    for k, v in kv:
        if k in kj:
            collected[k + ' [kJ/mol]'] = v
        elif k in jmol:
            collected[k + ' [J/(mol K)]'] = v
        else:
            collected[k] = v 

    pd.DataFrame(collected).to_excel('thermo.xlsx', index = False) 
