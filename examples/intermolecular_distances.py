#!/usr/bin/env python3

"""
File: intermolecular_distances.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Find intermolecular distances in ionic liquid-water mixtures.
Cation-anion, Cation-water and anion-water lengths are found. Special care
is taken to avoid alkyl interactions with neighbouring ions, but the C2 proton
of imidazolium cations is accounted for- in fact this will be the only atom
from imidazolium atoms contributing to these bonds.
"""

import chem_assistant as ca
import glob

files = glob.glob('*xyz') # or import the get_files() function from chem_assistant: get_files('.', 'xyz')


def imid_c2_h(atom):
    """
    Checks if a C2-H proton of imidazolium is found
    """
    connectors = {}
    if atom.symbol == 'H':        
        for alpha in atom.connected_atoms: # alpha = one atom away, beta = two away
            for beta in alpha.connected_atoms: 
                if beta.symbol not in connectors:
                    connectors[beta.symbol] = 1
                else:
                    connectors[beta.symbol] += 1
                if alpha.symbol == 'C' and 'N' in connectors and connectors['N'] == 2:
                    return True
    return False

def not_alkyl(atom):
    if atom.symbol == 'H':
        for a in atom.connected_atoms:
            if a.symbol == 'C':
                return False
    return True


def interatomic_dist(mol, i, j):
    """
    Returns interatomic distance between i and j and determines
    the molecules containing atoms i and j
    """
    i_mol = mol.fragments[i.mol]['name']
    j_mol = mol.fragments[j.mol]['name']
    dist = i.distance_to(j)
    if dist < 2.2 and not i.symbol == j.symbol == 'H': # two hydrogens were being found...
        if i_mol in ca.Molecule.Cations and j_mol is 'water':
            if not_alkyl(i) or imid_c2_h(i):
                # print('Cation-water', f'{i_mol} ({i.symbol}, {i.index})', f'water ({j.symbol}, {j.index})', f'{dist:.2f}')
                return 'Cation-Water', dist
        if i_mol in ca.Molecule.Anions and j_mol is 'water':
            if not_alkyl(i):
                # print('Anion-water', f'{i_mol} ({i.symbol}, {i.index})', f'water ({j.symbol}, {j.index})', f'{dist:.2f}')
                return 'Anion-Water', dist
        if i_mol in ca.Molecule.Cations and j_mol in ca.Molecule.Anions:
            if not_alkyl(i) or imid_c2_h(i):
                # print('Cation-Anion', f'{i_mol} ({i.symbol}, {i.index})', f'{j_mol} ({j.symbol}, {j.index})',f'{dist:.2f}')
                return 'Cation-Anion', dist
        if i_mol is 'water' and j_mol is 'water':
            return 'Water-Water', dist


bonds = {}
for f in sorted(files):
    atoms_already_considered = []
    bonds[f] = {}
    print(f)
    mol = ca.Molecule(using=f)
    mol.separate()
    for i in mol.coords:
        for j in mol.coords:
            if i.mol != j.mol:  # different fragments
                indices = sorted([i.index, j.index])
                if indices not in atoms_already_considered:
                    ret = interatomic_dist(mol, i, j)
                    if ret is not None:
                        bond, dist = ret
                        if bond not in bonds[f]:
                            bonds[f][bond] = [dist]
                        else:
                            bonds[f][bond].append(dist)
                        atoms_already_considered.append(indices)

    # if any are not found, still need to add an empty list
    for bond in ('Cation-Water', 'Anion-Water', 'Cation-Anion', 'Water-Water'):
        if bond not in bonds[f]:
            bonds[f][bond] = []
           
# Just pretty-printing in using the responsive_table() utility. 
# Could easily just export to a csv here with ca.write_csv_from_dict(bonds) or similar

files = []
cat_an = []
cat_w = []
an_w = []
w_w = []
bonding = {'Cation-Water': cat_w, 'Anion-Water': an_w, 'Cation-Anion': cat_an, 'Water-Water': w_w}

for file, val in bonds.items():
    num_rows_per_file = max(len(dists) for dists in val.values())
    for _ in range(num_rows_per_file):
        files.append(file)
    for b, dists in val.items():
        # make lists same length as file column
        to_add = num_rows_per_file - len(dists)
        for _ in range(to_add):
            dists.append('NA')
        for b_type, lst in bonding.items():
            if b == b_type:
                lst += dists

flattened = {'File': files, 'Cation-Anion': cat_an, 'Cation-Water': cat_w, 'Anion-Water': an_w, 'Water-Water': w_w}
ca.responsive_table(flattened, strings = [1])
ca.write_csv_from_dict(flattened, 'intermolecular_distances.csv')
