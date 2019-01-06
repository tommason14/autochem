#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['results_table', 'parse_results', 'thermochemistry', 'get_results_class',
'search_for_coords', 'get_h_bonds']

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
from ..core.molecule import Molecule
from ..core.thermo import thermo_data
from ..core.utils import (read_file, get_files, write_csv_from_dict, write_csv_from_nested)
from ..interfaces.gamess_results import GamessResults
from ..interfaces.psi_results import PsiResults
import os
import re
import csv

def search_for_coords(dir):
    """
    Recursively searched log/out files of optimisations for a successful equilibration- then writes to `equil.xyz`. If unsuccesful, writes to `rerun/rerun.xyz`, whilst also creating the corresponding input and job file.
    """
    
    for log in get_files(dir, ('.log', '.out')):
        r = get_results_class(log)
        if r.completed():
            if r.is_optimisation():
                print(f'{r.log}:\nFinding equilibrium coordinates...', end = " ")
                r.get_equil_coords()
                print()
        else:
            print(f"{log}: Not completed")


def get_results_class(log):
    """Return an instance of the desired class- |GamessResults|, |PsiResults|"""
    log_type = get_type(log)
    logs = {'gamess': GamessResults(log),
            'psi': PsiResults(log)}
    return logs.get(log_type, None)


def get_type(filepath):
    """
    Read in file, determine calculation type
    """
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

def parse_results(dir):
    """
    Used internally to parse log files for energies
    """
    output = []
    for log in get_files(dir, ('.out', '.log')):
        calc = get_results_class(log)
        try:
            if calc.completed(): #add provision for energies of opts only if equilibrium found
                if not calc.is_hessian():
                    print(calc.log)
                    data = calc.get_data()
                    output.append(data)
            else:
                print(f'{calc.log}: Incomplete')
        except AttributeError: # if log/out files are not logs of calculations
            continue
    return output    

def results_table(dir):
    """
    Prints energies of all log/out files in current and any sub directories to the screen, with the option of saving to csv.
    """
    output = parse_results(dir)
    print(f"{'File':^30s} | {'Path':^60s} | {'Basis':^8s} | {'HF':^15s} | {'MP2/SRS':^15s}")
    print('-'*140)
    for res in output:
        # need to convert empty HF to numeric
        f, p, b, hf, mp2 = res
        if hf == '':
            hf = 0
        print(f"{f:^30s} | {p:^60s} | {b:^8s} | {hf:^15.6f} | {mp2:^15.6f}")
    name = write_csv_from_nested(output, col_names = ('File', 'Path', 'Basis', 'HF', 'MP2/SRS'), return_name = True)
    return name # for use in other calculations (chem_assist -r uses this name)


def thermochemistry(dir):
    """
    Returns thermochemical data for all the relevant hessian log files in the given directory and
    subdirectories. Saves to csv file.
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
    write_csv_from_dict(collected)

def get_h_bonds(dir):
    """
    Searches the current and any subdirectories for xyz files, then attempts to split them into fragments, reporting any intermolecular bonding involving hydrogen-bonding atoms less than 2 Å apart. Prints to screen, with then option of saving to csv.

    TODO: Include the atoms bonded to the atom undergoing hydrogen-bonding. Two-fold benefit; can then disregard interactions involving alkyl chains, and can include angle information- internal angles of hydrogen bonds (connected-donor---acceptor) must be <= 45°
    """
    output = []
    for file in get_files(dir, ("xyz")):
        path, f = os.path.split(file)
        if not any((re.search('cation-?_?[0-9]*', f), re.search('anion-?_?[0-9]*', f), re.search('neutral-?_?[0-9]*', f))) and not any((f.split('_')[0] in Molecule.Anions, f.split('_')[0] in Molecule.Cations, f.split('_')[0] in Molecule.Neutrals)) and 'frags' not in path:
            # no frags in path of xyz- and no files named, cation_1.xyz, cation-1.xyz, anion, neutral
            # check for names in Anions, Cations, Neutrals- only want complexes
            print(file + '\n')
            mol = Molecule(using = file)
            mol.nfrags = int(input('Number of fragments: '))
            mol.separate()
            res = mol.find_h_bonds()
            for i in res:
                i.insert(0, f)
                i.insert(1, path)
            output += res
    write_csv_from_nested(output, col_names=('File', 'Path', 'Molecules', 'Length (Å)'))
