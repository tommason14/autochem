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
TODO:
- If geom_opt has completed, extract equilibrium coords: save in parent folder as equil.xyz & if no
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
import re
import csv


def get_results_class(log):
    """Return an instance of the desired class- |GamessResults|, |PsiResults|"""
    log_type = get_type(log)
    logs = {'gamess': GamessResults(log),
            'psi': PsiResults(log)}
    return logs.get(log_type, None)

def search_for_coords(dir):
    for log in get_files(dir, ('.log', '.out')):
        r = get_results_class(log)
        if r.completed():
            if r.is_optimisation():
                print(f'{r.log}: Finding equilibrium coordinates...')
                r.get_equil_coords()
        else:
            print(f"{log}: Not completed")

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
        # if r.get_fmo_level() != 0:
        #     energy = srs_output(r)
        # else:
        #     energy = (r.get_non_fmo(),) # Total Energy = .... 
        data = r.get_mp2_data()
    if type(r) == PsiResults:
        energy = srs_output(r)
    return data


def get_type(filepath):
    calc = ''
    with open(filepath, "r") as f:
        for line in f.readlines():
            if 'GAMESS' in line:
                calc = 'gamess'
                break
            elif 'PSI4' in line:
                calc = 'psi'
                break
    return calc

def gamess_run_type(filepath):
    fmo = False
    mp2 = False
    scs = False
    # fmo, mp2/srs, hf, dft?
    with open(filepath, "r") as f:
        for line in f.readlines():
            if 'FMO' in line:
                fmo = True
            elif 'MPLEVL' in line:
                mp2 = True
            elif 'SCS' in line:
                scs = True
            elif 'RUN TITLE' in line:
                break # by this point, all data required is specified
    return fmo, mp2, scs

def mp2_data(filepath, mp2_type):
    basis = ''
    HF = ''
    MP2 = ''
    with open(filepath, "r") as f:
        for line in f.readlines():
            if 'Euncorr HF' in line:
                HF = float(line.split()[-1])
            if 'INPUT CARD> $BASIS' in line:
                basis = line.split()[-2].split('=')[1]
            if f'E corr {mp2_type}' in line:
                MP2 = float(line.split()[-1])
    return basis, HF, MP2

def gamess_data(filepath):
    """Returns the last occurrence of FMO energies (FMO3 given if available), SRS and HF energies"""
    fmo, mp2, scs = gamess_run_type(filepath)
    if fmo and scs:
        basis, HF, MP2 = mp2_data(filepath, 'SCS')
    elif fmo and mp2 and not scs:
        basis, HF, MP2 = mp2_data(filepath, 'MP2')
    elif not fmo:
        val = ''
        basis = ''
        HF = ''
        MP2 = ''
        
        with open(filepath, "r") as f:
            for line in f.readlines():
                if 'INPUT CARD> $BASIS' in line:
                    basis = line.split()[-2].split('=')[1]
                if 'TOTAL ENERGY =' in line:
                    val = float(line.split()[-1])
        # What is total energy? SRS/MP2 or HF/DFT or something else?
        if scs or mp2:
            MP2 = val
        else:
            HF = val

    # more readable basis set
    change_basis = {'CCD'  : 'cc-pVDZ',
                    'CCT'  : 'cc-pVTZ',
                    'CCQ'  : 'cc-pVQZ',
                    'aCCD' : 'aug-cc-pVDZ',
                    'aCCT' : 'aug-cc-pVTZ',
                    'aCCQ' : 'aug-cc-pVQZ'}
    basis = change_basis.get(basis, basis) # default is the current value

    return basis, HF, MP2        


def psi_data(filepath):
    HF = ''
    opp = ''
    same = ''
    basis = ''
    MP2 = '' 
    with open(filepath, "r") as f:
        for line in f.readlines():
            if re.search('basis\s\w*(\-?\w*){1,2}$', line):
                basis = line.split()[-1]
            if 'Reference Energy          =' in line:
                HF = float(line.split('=')[1].split()[0].strip())
            elif 'Same-Spin Energy          =' in line:
                same = float(line.split('=')[1].split()[0].strip())
            elif 'Opposite-Spin Energy      =' in line:
                opp = float(line.split('=')[1].split()[0].strip())
    
    SRS = {}
    SRS['c_os'] = {}
    SRS['c_os']['ccd']         = 1.752
    SRS['c_os']['cc-pvdz']     = 1.752
    SRS['c_os']['cct']         = 1.64
    SRS['c_os']['cc-pvtz']     = 1.64
    SRS['c_os']['ccq']         = 1.689
    SRS['c_os']['cc-pvqz']     = 1.689
    SRS['c_os']['accd']        = 1.372
    SRS['c_os']['aug-cc-pvdz'] = 1.372
    SRS['c_os']['acct']        = 1.443
    SRS['c_os']['aug-cc-pvtz'] = 1.443
    SRS['c_os']['accq']        = 1.591
    SRS['c_os']['aug-cc-pvqz'] = 1.591

    SRS['c_ss'] = {}
    SRS['c_ss']['ccd']         = 0
    SRS['c_ss']['cc-pvdz']     = 0
    SRS['c_ss']['cct']         = 0
    SRS['c_ss']['cc-pvtz']     = 0
    SRS['c_ss']['ccq']         = 0
    SRS['c_ss']['cc-pvqz']     = 0
    SRS['c_ss']['accd']        = 0
    SRS['c_ss']['aug-cc-pvdz'] = 0
    SRS['c_ss']['acct']        = 0
    SRS['c_ss']['aug-cc-pvtz'] = 0
    SRS['c_ss']['accq']        = 0
    SRS['c_ss']['aug-cc-pvqz'] = 0


    c_os = SRS['c_os'][basis.lower()]
    c_ss = SRS['c_ss'][basis.lower()]
    
    MP2 = HF + c_os * opp + c_ss * same
    
    return basis, HF, MP2

    

def get_data(log):
    path, file = os.path.split(log)
    calctype = get_type(log)
    if calctype == 'gamess':
        basis, HF, MP2 = gamess_data(log)
    elif calctype == 'psi':
        basis, HF, MP2 = psi_data(log)
    return file, path, basis, HF, MP2


def parse_results(dir):
    output = []
    cwd = os.getcwd()
    for log in get_files(dir, ('.out', '.log')):
        r = get_results_class(log)
        if r.completed(): #add provision for energies of opts only if complete
            if not r.is_hessian():
                print(log)
                data = get_data(log)
                output.append(data)
    return output    
    #     r = get_results_class(log)
    #     filetype = get_type(log)
    #     try:
    #         if r.completed():
    #             basis = r.get_basis()
    #             if r.is_optimisation() and r.opt_found_equil_coords() or r.is_spec(): # long time!!
    #                 print(f'{r.log}: Finding energies')
    #                 f = fetch_energies(r)
    #                 print(f)
    #                 # adding to csv/df/database
    #                 # output.append((log, filetype, basis) + f)
    #             elif r.is_hessian():
    #                 print(f'{r.log}: Hessian calc')
    #         else:
    #             # incomplete cases
    #             r.get_error()
    #     except AttributeError:
    #         continue
    # # return output

def write_csv(data):
    done = False
    while not done:
        to_file = input('Print to csv? [y/n] ')
        if to_file.lower() in ('y', 'n'):
            done = True
            if to_file.lower() == 'y':
                filename = input('Filename: ')
                with open(filename, "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(('File', 'Path', 'Basis', 'HF', 'MP2/SRS'))       
                    writer.writerows(data)
        else:   
            print("Please select 'y' or 'n'")

def results_table(dir):
    output = parse_results(dir)
    print(f"{'File':^30s} | {'Path':^60s} | {'Basis':^8s} | {'HF':^15s} | {'MP2/SRS':^15s}")
    print('-'*140)
    for res in output:
        # need to convert to floats for printing 
        f, p, b, hf, mp2 = res
        if hf == '':
            hf = 0
        print(f"{f:^30s} | {p:^60s} | {b:^8s} | {hf:^15.6f} | {mp2:^15.6f}")#.format(f, p, b, hf, mp2))
    write_csv(output) 


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
        r = get_results_class(log)
        if r.completed():
            if r.is_hessian():
                print(f'Thermo data for {log}')
                res = thermo_data(r.log) # run fortran script
                res['File'] = log
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
    write_thermo(collected)

def write_thermo(data):
    """Write to file from dictionary"""
    done = False
    while not done:
        to_file = input('Print to csv? [y/n] ')
        if to_file.lower() in ('y', 'n'):
            done = True
            if to_file.lower() == 'y':
                filename = input('Filename: ')
                with open(filename, "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(data.keys())       
                    writer.writerows(zip(*[data[key] for key in data.keys()]))
        else:   
            print("Please select 'y' or 'n'")
