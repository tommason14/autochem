#!/usr/bin/env python3
from chem_assistant import read_file, get_files, write_csv_from_nested, Molecule, Atom
import re
import os

cm5_results = []

files = get_files('.', 'cm5.out')
regex = "\s*[0-9]{1,3}\s*[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){6}"
inp_regex = '^[A-Z][a-z]?(\s*-?[0-9]*.[0-9]*){3}'

found = False
for file in files:
    path, f = os.path.split(file)
    print(path)
    inputfile = file.replace('out', 'inp')
    data = []
    coords = []
    parsed = []
    lines = read_file(file)
    for line in lines:
        if 'Q-H        S-H        Dx         Dy         Dz        Q-CM5' in line:
            found = True
        elif 'Tot' in line:
            found = False
        if found:
            # index, sym, hirshfeld charge, spin density, dipole_x, dipole_y, dipole_z, cm5 charge
            if re.search(regex, line):
                res = line.split()
                data.append([path] + res[:3] + res[4:])
    
    # need coords
    for line in read_file(inputfile):
        if re.search(inp_regex, line):
            sym, x, y, z = line.split()
            coords.append(Atom(symbol = sym, coords = (x, y, z)))

    mol = Molecule(atoms = coords)
    mol.separate()
    
    # # naming- i.e. c4mim_1, c4mim_2, water_1
    # print(mol.fragments)
    #
    # mols = {}    
    #
    # for data in mol.fragments.values():
    #     if data['name'] not in mols:
    #         mols[data['name']] = 1
    #     else:
    #         mols[data['name']] += 1
    # print(mols)
    #
    # molecules = []
    # for mole, num in mols.items():
    #     for i in range(1, num + 1):
    #         molecules.append(mole + '_' + str(i))
    # print(molecules)

    for atom in mol.coords:
        # print(atom.mol, molecules[atom.mol], atom.index)
        parsed.append([atom.x, atom.y, atom.z, mol.fragments[atom.mol]['name'] + '_' + str(atom.mol)])
    

    for line in zip(data, parsed):
        cm5_results.append((*line[0], *line[1]))

write_csv_from_nested(cm5_results, col_names = ('Path', 'Index', 'Element', 'Hirshfeld', 'Dx', 'Dy',
'Dz', 'CM5', 'Rx', 'Ry', 'Rz', 'Fragment'))
