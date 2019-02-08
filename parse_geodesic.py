#!/usr/bin/env python3

from chem_assistant import (Atom, Molecule,
read_file, get_files, write_csv_from_nested)
import os
import sys
import re

atom_regex = '^\s[A-Za-z]{1,2}\s*[0-9]*.[0-9]*(\s*-?[0-9]*.[0-9]*){3}$'
charge_regex = '^\s[A-Za-z]{1,2}(\s*-?[0-9]*.[0-9]*){2}$'

results = []

files = get_files('.', 'log')

for logfile in files:
    path, filename = os.path.split(logfile)
    print(path)
    inpfile = logfile[:-3] + 'inp'

    res = []
    assigned = []

    for line in read_file(inpfile):
        if re.search(atom_regex, line):
            sym, atnum, x, y, z = line.split()
            x, y, z = map(float, (x, y, z))
            res.append([path, Atom(sym, coords = (x, y, z))]) # new key for each coord
    found = False
    counter = 0
    for line in read_file(logfile):
        if 'NET CHARGES:' in line:
            found = True
        if 'RMS DEVIATION' in line:
            break
        if found: 
            if re.search(charge_regex, line):
                res[counter].append(float(line.split()[1]))
                counter += 1

    coords = [atom[1] for atom in res]
    mol = Molecule(atoms = coords)
    mol.separate()
    for i, atom in enumerate(mol.coords):
        path, _, geodesic_charge = res[i]
        index = i + 1
        results.append([path, index, atom.symbol, geodesic_charge, atom.x, atom.y, atom.z,
f"{mol.fragments[atom.mol]['name']}_{atom.mol}"])

write_csv_from_nested(results, col_names = ('Path', 'Index', 'Element', 'Geodesic', 'Rx', 'Ry', 'Rz', 'Fragment'))
