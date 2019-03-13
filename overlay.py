#!/usr/bin/env python3

from chem_assistant import Molecule, Atom, write_xyz
import sys

if len(sys.argv) != 3:
    err = """
Error: Incorrect argument calls
Syntax: overlay.py mol1.xyz mol2.xyz
"""
    sys.exit(err)

file1, file2 = sys.argv[1:]

mol1 = Molecule(file1)
mol2 = Molecule(file2)


# define a reference point


# make the first atom have position (0,0,0)
# with all other atoms relative to that
# find vector from atom1 to (0,0,0), then move all atoms by that vector

vector_to_move_by = mol1.coords[0].vector_to((0,0,0))

for index, atom in enumerate(mol1.coords):
    new_coords = tuple(i + j for i,j in zip(atom.coords,vector_to_move_by))
    mol1.coords[index] = Atom(atom.symbol, coords = new_coords)

write_xyz(mol1.coords, f"{file1.rsplit('.')[0]}_moved.xyz")

vector_to_move_by2 = mol2.coords[0].vector_to((0,0,0))

for index, atom in enumerate(mol2.coords):
    new_coords = tuple(i + j for i,j in zip(atom.coords,vector_to_move_by2))
    mol2.coords[index] = Atom(atom.symbol, coords = new_coords)

write_xyz(mol2.coords, f"{file2.rsplit('.')[0]}_moved.xyz")
