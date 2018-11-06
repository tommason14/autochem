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
from meta import Meta
import numpy as np
import math
from errors import *
import itertools

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

    def __init__(self, using = None):
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
        """Returns molecular mass in g mol⁻¹"""
        formula = self.formula(as_dict=True)
        mass = 0
        for element, number in formula.items():
            mass += PT.get_mass(element) * number
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


    def split(self, cutoff):
        """Assigns each atom in ``self.coords`` to a different fragment- very useful for FMO calculations."""

        # Loop over each atom in turn, and compare with every other atom in the list. If the distance between atoms is less than a cutoff distance, then consider them to be in the same fragment. Consider bonds as appropriate.
        
        mol = 0
        for i, atom_i in enumerate(self.coords):
            connected = False
            for j, atom_j in enumerate(self.coords):
                if i != j:
                    dist = atom_i.distance_to(atom_j)
                    if dist < cutoff:
                        connected = True
                        #connected
                        #same molecule 
                        if atom_i.mol is None and atom_j.mol is None:
                        #if not already assigned:
                            atom_i.connected_atoms.append(atom_j)
                            atom_j.connected_atoms.append(atom_i)
                            atom_i.bonds.append(Bond(atom_i, atom_j))
                            atom_j.bonds.append(Bond(atom_j, atom_i))
                            atom_i.mol, atom_j.mol = mol, mol
                            mol += 1
                        elif atom_i.mol is not None and atom_j.mol is None:
                            # only assign to j
                            # need to assign the number of i to j, not a new value
                            atom_j.mol = atom_i.mol
                            atom_j.connected_atoms.append(atom_i)
                            atom_j.bonds.append(Bond(atom_j, atom_i))
                        elif atom_j.mol is not None and atom_i.mol is None:
                            # only assign to i
                            atom_i.mol = atom_j.mol
                            atom_i.connected_atoms.append(atom_j)
                            atom_i.bonds.append(Bond(atom_i, atom_j))
                        elif atom_i.mol is not None and atom_j.mol is not None:
                            """MAY NEED TO ADD SOME BONDS HERE LATER
                            FOR THE HYDROGEN BOND DISTANCES ETC..."""
                            # if within cutoff, but assigned to different fragments, need to re-assign all atoms of one number to another
                            if atom_i.mol != atom_j.mol:
                                #look up a number to replace
                                num = atom_j.mol # i.e. get the 4s
                                for atom in self.coords:
                                    if atom.mol == num: # when atom.mol is 4
                                        atom.mol = atom_i.mol # reassign to 3
            if not connected:
                # assign its own fragment- say a halide anion
                atom_i.mol = mol
                mol +=1
        
    def collect_frags(self):
        """Collating all atoms of fragments into a list, assigned to a key of a dictionary. 
        Moving from a 1-dimensional list of coordinates to a dictionary of fragments"""

        self.frags = {}
        for atom in self.coords:
            if atom.mol not in self.frags:
                self.frags[atom.mol] = [atom]
            else:
                self.frags[atom.mol].append(atom)

    def check_db(self):
        """Checks fragments for a match in the database"""

        symbols = {} # not tied to an instance, only required when checking the database
        for frag, atoms in self.frags.items():
            symbols[frag] = [atom.symbol for atom in atoms]
        
        def check_dict(molecules_dict, symbols_dict, db):
            """Checking molecule database for a match, and returning the required attributes- name, atoms in molecule, type of molecule, charge, multiplicity"""
            if db == Molecule.Anions:
                charge = -1
                mult = 1
                mol_type = 'anion'
            elif db == Molecule.Cations:
                charge = 1
                mult = 1
                mol_type = 'cation'
            elif db == Molecule.Neutrals:
                charge = 0
                mult = 1
                mol_type = 'neutral'

            for name, atom_list in db.items():
                for sym, molecule in symbols_dict.items():
                    if sorted(molecule) == sorted(atom_list): #sorts in place, but no overwriting of variable
                        data = {
                                "type": mol_type,
                                "name" : name,
                                "atoms": self.frags[sym], #not symbols, but the atom objects
                                "charge": charge,
                                "multiplicity": mult
                            }
                        molecules_dict[sym] = data
            return molecules_dict

        self.fragments = {}
        for db in (Molecule.Anions, Molecule.Cations, Molecule.Neutrals):
            self.fragments = check_dict(self.fragments, symbols, db)
    
    def print_frags(self):
        for frag, data in self.fragments.items():
            print(f"{data['type'].capitalize()} found: {data['name']}")

    def remove_assignments(self):
        for atom in self.coords:
            atom.mol = None
    
    def separate(self):
        """Separates coordinates into specific fragments using a cutoff distance. Note this function only works with intermolecular fragments and cannot split molecules on bonds."""
        cutoff = 1.7
        correct = False
        while not correct:
            self.split(cutoff)
            self.collect_frags()
            self.check_db()
            self.print_frags()
            check = input(f'Does your system contain {len(self.fragments)} fragments? [y/n] ')
            if check.lower() == 'y':
                correct = True
            if check.lower() == 'n':
                self.remove_assignments()
                cutoff = float(input('Cutoff distance (Å): '))
            if check.lower() not in ('y', 'n'):
                print("Please choose 'y' or 'n'")

    #format for fmo 
    def gamess_format(self):
        """Returns a tuple of strings, for the INDAT and ICHARG blocks of GAMESS FMO calculations"""
        index_dict = {}
        for index, atom in enumerate(self.coords):
            if atom.mol not in index_dict:
                index_dict[atom.mol] = [index + 1]
            else:
                index_dict[atom.mol].append(index + 1)
        
        indat = ""
        for frag, lst in index_dict.items():
            indat +=  f"0,{lst[0]},{lst[-1]},\n"
        indat += "0"

        # as python dicts are arranged in numerical / alphabetical order, the order of mol in self.coords and the order of fragments in self.frags are the same
        icharg = ""
        for frag, data in self.fragments.items():
            icharg += f"{data['charge']},"
        icharg = icharg[:-1] #rm trailing comma
        return (indat, icharg)



