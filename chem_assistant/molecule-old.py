#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File: molecule.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Holds class representing a molecule
"""

from periodic_table import PeriodicTable as PT
from atom import Atom
from bond import Bond
import numpy as np
import math
from errors import *


class Molecule:
    """Class implementing the concept of a molecular system. Crucial is the separation of molecules in a given system, with output for Fragment Molecular Orbital calculations"""

    Anions = {"br" : ["Br"]}
    Anions["cl"] = ["Cl"]
    Anions["bf4"] = ['B', 'F', 'F', 'F', 'F']
    Anions["dca"] = ['N', 'C', 'N', 'C', 'N']
    Anions["pf6"] = ['F', 'P', 'F', 'F', 'F', 'F', 'F']
    Anions["mes"] = ['S', 'O', 'O', 'O', 'C', 'H', 'H', 'H']
    Anions["ntf2"] = ['F', 'F', 'F', 'F', 'F', 'N', 'S', 'S',
                        'O', 'O', 'O', 'O', 'C', 'C', 'F']
    Anions["tos"] = ['C', 'C', 'C', 'C', 'H', 'H', 'H', 'H',
                        'H', 'H', 'H', 'S', 'O', 'O', 'O', 'C', 'C', 'C']
    Anions["h2po4"] = ['H', 'H', 'P', 'O', 'O', 'O', 'O']
    Anions["acetate"] = ['C','H', 'H', 'H', 'C', 'O', 'O']

    Cations = {"c1mim": ['C', 'N', 'C', 'N',
                        'C', 'C', 'C', 'H',
                        'H', 'H', 'H', 'H',
                        'H', 'H', 'H', 'H']}
    Cations["c1mpyr"] = ['C', 'C', 'C', 'N', 'C', 'C', 'C',
                            'H', 'H', 'H', 'H', 'H', 'H', 'H',
                            'H', 'H', 'H', 'H', 'H', 'H', 'H']

    Cations["c2mim"] = ['C', 'N', 'C', 'N', 'C',
                        'C', 'C', 'C', 'H', 'H',
                        'H', 'H', 'H', 'H', 'H',
                        'H', 'H', 'H', 'H']
    Cations["c2mpyr"] = ['N', 'C', 'C', 'C', 'C', 'C',
                        'C', 'C', 'H', 'H', 'H', 'H', 'H',
                        'H', 'H', 'H', 'H', 'H', 'H', 'H',
                        'H', 'H', 'H', 'H']

    Cations["c3mim"] = ['N', 'C', 'N', 'C', 'C', 'C', 'C', 'H',
                            'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H',
                            'H', 'H', 'H', 'H', 'H', 'H']
    Cations["c3mpyr"] = ['N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
                            'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
                            'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']

    Cations["c4mim"] = ['C', 'N', 'C', 'C', 'N', 'C', 'C', 'H', 'C', 'C',
                            'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
                            'H', 'H', 'H', 'H', 'H', 'H']
    Cations["c4mpyr"] = ['N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
                            'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
                            'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
                            'H', 'H', 'H']
    Cations['choline'] = ['N', 'C', 'H', 'H', 'H', 'C', 'H', 'H',
                        'H', 'C', 'H', 'H', 'H', 'C', 'H', 'H',
                        'C', 'H', 'H', 'O', 'H']

    Neutrals = {"nh3" : ['N', 'H', 'H', 'H']}
    Neutrals['water'] = ['H', 'H', 'O']

    def __init__(self, using):
        if using is not None:
            self.coords = self.read_xyz(using)
            #list of Atom objects, more useful than list of coordinates
        self.bonds = []

    def __repr__(self):
        return f'Molecule of {len(self.coords)} atoms'

    def __iter__(self):
        return iter(self.coords)

    def formula(self, as_dict = False, as_latex = False, as_html = False):
        """Returns the molecular format in a variety of formats using keyword arguments:
        * *as_dict* -- returns a python dictionary
        * *as_latex* -- returns a latex expression using \textsubscript{} notation (equivalent to $_{value}$, but the math expression renders in a different font)
        * *as_html* -- returns a html expression using <sub> tags

        If no keyword arguments are passed, a regular string is returned with no formatting i.e. C8H18

        """

        formula  = {}
        for atom in self.coords:
            if atom.symbol not in formula:
                formula[atom.symbol] = 1
            else:
                formula[atom.symbol] += 1
        if as_dict:
            return formula
        if as_latex:
            string = ""
            for sym, number in sorted(formula.items()):
                string += f"{sym}\\textsubscript{{{number}}}"
            return string
        if as_html:
            string = ""
            for sym, number in sorted(formula.items()):
                string += f"{sym}<sub>{number}</sub>"
            return string
        else:
            string = ""
            for sym, number in sorted(formula.items()):
                string += f"{sym}{number}"
            return string

    def mass(self):
        formula = self.formula(as_dict=True)
        mass = 0
        for element, number in formula.items():
            mass += PT.get_mass(element) * number
        # return f"{mass:.2f} g mol \u207b\u00b9"
        return f"{mass:.2f} g mol⁻¹"

    def read_xyz(self, using):
        """Reads coordinates of an xyz file and return a list of |Atom| objects, one for each atom"""
        coords = []
        with open(using, "r") as f:
            for coord in f.readlines()[2:]:
                line = coord.split()
                for val in PT.ptable.values():
                    if line[0] == val[0]:
                        coords.append(Atom(line[0], coords = tuple(float(i) for i in line[1:4])))
        return coords

    def write_xyz(self, atoms, filename = None):
        """Writes an xyz file using a list of |Atom| instances"""
        if filename is None:
            raise ValueError('write_xyz: Must give a path to the output file')
        else:
            with open(filename, "w") as f:
                f.write(str(len(atoms)) + '\n\n')
                for atom in atoms:
                    f.write(f"{atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f} \n")


    def separate(self):
        """Assigns each atom in ``self.coords`` to a different fragment- very useful for FMO calculations"""

        def atom_distances(self):
            """Returns numpy array of interatomic distances between all atoms of the molecule"""
            dimension = len(self.coords)
            dist = np.zeros((dimension, dimension))
            for i, atom_i in enumerate(self.coords):
                for j, atom_j in enumerate(self.coords):
                    if atom_i != atom_j:
                        dist[i, j] = atom_i.distance_to(atom_j)
            return dist

        def assign_group(self, distances, cutoff = None):
            """Assigns a fragment to all atoms in the molecule. Updates the ``mol`` attribute of |Atom| and |Bond| instances with a number, and different numbers represent different fragments"""
            if cutoff == None:
                cutoff = 1.7
            num = 0
            for i, atom_i in enumerate(self.coords):
                con = False
                for j, atom_j in enumerate(self.coords):
                    if atom_i != atom_j:
                        if distances[i, j] < cutoff:
                            con = True
                            if distances[i, j] < (PT.get_radius(atom_j.symbol) + PT.get_radius(atom_j.symbol)): 
                                # assign same number to each atom, and to bonds between
                                # those atoms
                                
                                if atom_j not in atom_i.connected_atoms:
                                    atom_i.connected_atoms.append(atom_j)
                                    atom_j.connected_atoms.append(atom_i)
                                    atom_i.bonds.append(Bond(atom_i, atom_j, mol = num))
                                    atom_j.bonds.append(Bond(atom_j, atom_i, mol = num))
                                if atom_i.mol is None and atom_j.mol is None:
                                    atom_i.mol, atom_j.mol = num, num
                                    num += 1
                                elif atom_i.mol is not None and atom_j.mol is not None:
                                    if atom_i.mol != atom_j.mol:
                                        # atoms assigned to different groups, but they are in the same fragment-- correction for this
                                        lookup = atom_j.mol
                                        for atom in self.coords:
                                            if atom.mol is lookup:
                                                atom.mol = atom_i.mol
                                        # look for all of one fragment, and assign them all to the other number
                                        # i.e. 1,2,1,1,1 -> look for 2, assign as 1 -> 1,1,1,1,1

                                # if one atom has no assigned number, but is within cutoff
                                elif atom_i.mol is not None and atom_j.mol is None:
                                    atom_j.mol = atom_i.mol
                                elif atom_j.mol is not None and atom_i.mol is None:
                                    atom_i.mol = atom_j.mol
                if not con:
                    atom_i.mol = num
                    num += 1

        def grouping(self):
            """Returns a dictionary of fragments with a list of atoms as their value"""
            fragments = {}
            for atom in self.coords:
                if atom.mol not in fragments:
                    fragments[atom.mol] = [atom]
                else:
                    fragments[atom.mol].append(atom)
            return fragments

        def check_db(self, frags):
            """Returns a dictionary of fragments, complete with names, multiplicities and charges, ready for output either to individual files, or as input for a GAMESS FMO calculation"""
            molecules = {}

            symbols = {}
            for frag in frags:
                for atom in frags[frag]:
                    if frag not in symbols:
                        symbols[frag] = [atom.symbol]
                    else:
                        symbols[frag].append(atom.symbol)

            def check_dict(molecules_dict, symbols_dict, frags, db, charge, mult, mol_type):
                """Checking molecule database for a match, and returning the required attributes- name, atoms in molecule, type of molecule, charge, multiplicity"""
                for name, atom_list in db.items():
                    for sym, molecule in symbols_dict.items():
                        if sorted(molecule) == sorted(atom_list): #sorts in place, but no overwriting of variable
                            molecules_dict[sym] = {
                                "type": mol_type,
                                "name" : name,
                                "atoms": frags[sym], #not symbols, but the atom objects
                                "charge": charge,
                                "multiplicity": mult
                            }
                return molecules

            molecules = check_dict(molecules, symbols, frags, Molecule.Anions, -1, 1, 'anion')
            molecules = check_dict(molecules, symbols, frags, Molecule.Cations, 1, 1, 'cation')
            molecules = check_dict(molecules, symbols, frags, Molecule.Neutrals, 0, 1, 'neutral')

            return molecules



        def assign_frags(distances, cutoff=None):
            assignments = assign_group(self, distances, cutoff)
            groups = grouping(self)
            frags = check_db(self, groups)
            return frags

        def print_fragments(frags):
            for frag, meta in frags.items():
                print(f"{meta['type'].capitalize()} found: {meta['name']}")
            try:
                call = input(f"Does your system have {len(frags)} fragments? [Y/N] ").lower()
                if call not in ('y', 'n'):
                    print("Please choose 'Y' or 'N'")
                else:
                    if call == 'y':
                        correct = True
                    if call == 'n':
                        correct = False
            except ValueError:
                print("Please choose 'Y' or 'N'")
                print_fragments(frags)
            return correct


        def separate_main():
            distances = atom_distances(self)
            frags = assign_frags(distances)
            correct = print_fragments(frags)
            while not correct:
                no_value = True
                while no_value:
                    try:
                        cutoff = float(input('Cutoff distance (Å): '))
                        no_value = False
                    except ValueError:
                        print('Please enter a numerical value')
                frags = assign_frags(distances, cutoff)
                correct = print_fragments(frags)

        separate_main()
