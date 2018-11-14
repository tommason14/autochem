#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File: utils.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Functions useful in multiple scenarios 
"""
__all__ = ['read_file', 'get_type']

def read_file(file):
    with open(file, "r") as f:
        for line in f.read():
            yield line

def get_type(file):
    for line in read_file(file):
        if 'PSI4' in line:
            return 'psi4'
        elif 'GAMESS' in line:
            return 'gamess'
        # extend to lammps

def write_xyz(atoms, filename = None, defined_atoms = False):
    """Writes an xyz file using a list of |Atom| instances"""
    if filename is None:
        raise ValueError('write_xyz: Must give a path to the output file')
    else:
        with open(filename, "w") as f:
            f.write(str(len(atoms)) + '\n\n')
            for atom in atoms:
                if not defined_atoms:
                    f.write(atom + '\n')
                else:
                    f.write(f"{atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f} \n")
