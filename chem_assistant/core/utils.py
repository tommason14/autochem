#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .atom import Atom

__all__ = ['read_file', 'get_type']

def read_file(file):
    with open(file, "r") as f:
        for line in f.readlines():
            yield line

def get_type(file):
    for line in read_file(file):
        if 'PSI4' in line:
            return 'psi4'
        elif 'GAMESS' in line:
            return 'gamess'
        # extend to lammps

def write_xyz(atoms, filename = None):
    """Writes an xyz file using a list of |Atom| instances, or just a list of regular coordinates,
with or without atomic numbers."""
    if filename is None:
        raise ValueError('write_xyz: Must give a path to the output file')
    else:
        with open(filename, "w") as file:
            file.write(str(len(atoms)) + '\n\n')
            for atom in atoms:
                if type(atom) is not Atom:
                    parts = atom.split()
                    if len(parts) > 4: # includes atomic nums
                        sym, _, x, y, z = parts
                    else:
                        sym, x, y, z = parts
                    x, y, z = float(x), float(y), float(z)
                    file.write(f"{sym:5s} {x:>13.10f} {y:13.10f} {z:13.10f} \n")
                else:
                    file.write(f"{atom.symbol:5s} {atom.x:>13.10f} {atom.y:>13.10f} {atom.z:>13.10f} \n")
