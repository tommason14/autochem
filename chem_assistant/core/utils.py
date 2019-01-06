#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv

from .atom import Atom
from .periodic_table import PeriodicTable as PT

__all__ = ['read_file', 'get_type', 'write_xyz', 'get_files', 'module_exists', 'sort_elements', 'write_csv_from_dict', 'write_csv_from_nested']

def read_file(file):
    with open(file, "r") as f:
        try:
            for line in f.readlines():
                yield line
        except UnicodeDecodeError:
            pass

def get_type(file):
    for line in read_file(file):
        if 'PSI4' in line:
            return 'psi4'
        elif 'GAMESS' in line:
            return 'gamess'
        # extend to lammps

def write_xyz(atoms, filename = None):
    """Writes an xyz file using a list of |Atom| instances, or just a list of regular coordinates,
with or without atomic numbers.
"""
    if filename is None:
        raise ValueError('write_xyz: Must give a path to the output file')
    else:
        with open(filename, "w") as file:
            file.write(str(len(atoms)) + '\n\n')
            for atom in atoms:
                if type(atom) is not Atom:
                    parts = atom.split()
                    if len(parts) > 4: # includes atomic nums
                        sym, *_, x, y, z = parts
                    else:
                        sym, x, y, z = parts
                    x, y, z = float(x), float(y), float(z)
                    file.write(f"{sym:5s} {x:>15.10f} {y:15.10f} {z:15.10f} \n")
                else:
                    file.write(f"{atom.symbol:5s} {atom.x:>15.10f} {atom.y:>15.10f} {atom.z:>15.10f} \n")


def get_files(directory, ext):
    """Accepts a tuple of file extensions, searches in all subdirectories of the directory given for relevant files. Returns a list of files with their relative path to the directory passed in.

    Usage:
        >>> for filepath in get_files('.', ("log", "out")):
        >>>     parse_file(filepath)
    """
    fileLst = []
    for path, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(ext) and file != 'freq.out': # freq.out used for thermo calculations with the fortran code
                fileLst.append(os.path.join(path, file))
    return fileLst

def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True

def sort_elements(lst):
    """
    Sort a list of |Atom| objects by atomic number. 
    Returns a list of tuples- [(symbol, atomic number), (symbol, atomic number), ...]

    TODO: Extend to giving back the objects- more useful than just for formatting of symbols
    """
    els = []
    elements = set([atom.symbol for atom in lst])
    for i in elements:
        els.append((i, float(PT.get_atnum(i))))
    sorted_els = sorted(els, key = lambda val: val[1])
    return sorted_els

def write_csv_from_dict(data):
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
                    content = zip(*[data[key] for key in data.keys()])
                    writer.writerows(content) 
        else:   
            print("Please select 'y' or 'n'")

def write_csv_from_nested(data,*,col_names = None, return_name = False):
    """
    Write to csv from nested data structure; list of tuples, list of lists. 
    
    NB: requires a list or tuple of column names passed to the `col_names` parameter
    """
    if type(col_names) not in (list, tuple):
        raise AttributeError('Must pass in column names as a list or tuple of values')

    done = False
    while not done:
        to_file = input('Print to csv? [y/n] ')
        if to_file.lower() in ('y', 'n'):
            done = True
            if to_file.lower() == 'y':
                filename = input('Filename: ')
                with open(filename, "w", encoding = 'utf-8-sig') as f:
                    writer = csv.writer(f)
                    writer.writerow(col_names)       
                    writer.writerows(data)
                if return_name:
                    return filename
        else:   
            print("Please select 'y' or 'n'")

def search_dict_recursively(d):
    ret = {}
    for k, v in d.items():
        if isinstance(v, dict):
            ret[k] = search_dict_recursively(v)
        else:
            ret[k] = v
    return ret
