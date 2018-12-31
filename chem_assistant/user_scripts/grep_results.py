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
    subdirectories.
    Usage:
        >>> thermochemistry('.')
    Returns:
        Thermodynamics for anion_0.out

        Thermochemistry was done at      298.15 (k) and 1(Atm)

        Total vibr frequencies                                           = 15
        Imaginary frequencies                                           =  0
        Low freqs (<300 cm^-1) with                      special treatment  =  0

        Scaling factors
        ZPVE =   1.0000
        TC   =   1.0000
        Svib =   1.0000

           Frequency     ZPVE(vi)      H(vi)      S(vi)
          170.5510         1.020109         1.597206        10.164604    HO
          266.5580         1.594352         1.217327         6.771481    HO
          333.6590         1.995700         0.996996         5.197825    HO
          401.0410         2.398729         0.809558         4.011760    HO
          430.5090         2.574985         0.737355         3.585659    HO
          479.1500         2.865919         0.630098         2.980525    HO
          490.0140         2.930899         0.608059         2.860051    HO
          746.1180         4.462723         0.250597         1.070728    HO
          816.4150         4.883188         0.193758         0.813204    HO
         1024.9390         6.130423         0.087821         0.353895    HO
         1045.1020         6.251023         0.081192         0.326142    HO
         1090.3080         6.521412         0.068016         0.271373    HO
         1272.6760         7.612202         0.032829         0.128019    HO
         3689.2490        22.066345         0.000001         0.000003    HO
         3723.2310        22.269600         0.000001         0.000002    HO

        ZPVE                        =            95.57761 kJ/mol

        Thermal correction in kJ/mol

        TC harmonic oscillator      =            17.22658

             Entropies in J/(mol K)
        S electronic                =             0.00000
        S translational             =           161.18860
        S rotational                =           102.72824
        S vibrational HO            =            38.53527

        Stotal Harmonic oscillator  =           302.45211

        TC-TdeltaS in kJ/mol        =           -72.94952
    """
    for log in get_files(dir, ('.log', '.out')):
        _, r = get_results_class(log)
        if r.completed():
            if r.is_hessian():
                thermo_data(r.log)
