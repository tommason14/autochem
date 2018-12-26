#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .periodic_table import PeriodicTable as PT
from .atom import Atom
from .bond import Bond
from .meta import Meta
from .errors import *

import numpy as np
import math
import itertools

__all__ = ['Molecule']

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


    def __init__(self, using = None, atoms = None, nfrags = None):
        if using is not None:
            self.coords = self.read_xyz(using)
            #list of Atom objects, more useful than list of coordinates
            for index, atom in enumerate(self.coords):
                atom.index = index + 1 #so first atom has index 1
        self.bonds = []
        if atoms is not None and using is None:
            self.coords = atoms

        self.nfrags = nfrags

        if hasattr(self, 'coords'):
        # self.complex used in input files
        # assuming a neutral closed shell system
            
            self.complex = {
                "type": "complex",
                "name": "complex",
                "atoms": self.coords,
                "charge": 0,
                "mult": 1,
                "elements": sort_elements(self.coords)
            }

    def __repr__(self):
        els = [i[0] for i in self.complex['elements']]
        if not hasattr(self, 'fragments'):
            return f'Molecule of {len(self.coords)} atoms. Elements: {els}'
        return self._repr()
    
    def _repr(self): 
        """Repr for system after fragmentation"""
        frags = {}
        for frag in self.fragments.values():
            if frag['name'] not in frags:
                frags[frag['name']] = 1
            else:
                frags[frag['name']] += 1
        string = f"Molecule of {len(self.fragments)} fragments, {len(self.coords)} atoms.\nFragments:\n"
        for frag, num in frags.items():
            string += f'    {num} x {frag}\n'
        return string[:-1] #lose last newline char

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
                string += f"{sym}\textsubscript{{{number}}}"
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


    def split(self):
        """Assigns each atom in ``self.coords`` to a different fragment- very useful for FMO calculations.
        Note: Only works if all intramolecular bonds are shorter than all intermolecular bonds."""

        self.mol_dict = {}
        for atom in self.coords:
            self.mol_dict[atom.index] = Molecule(atoms = [atom]) # 1-atom molecules

        def min_dist(mol_i, mol_j):
            """Finds the minimum distance between atoms in two different molecules"""
            distances = []
            for atom_i in mol_i:
                for atom_j in mol_j:
                    distances.append(atom_i.distance_to(atom_j))
            return min(distances)

        def distances():
            """Returns a list of lists of [i, j, min_dist(i, j)] for each combination of atoms in
the system"""
            dist_list = []
            # for every combination of atoms in the list, find their distance
            for i, j in itertools.combinations(self.mol_dict.keys(), 2):
                minimum_dist = min_dist(self.mol_dict[i], self.mol_dict[j])
                dist_list.append([i, j, minimum_dist])
            return dist_list

        def collect_atoms():
            # modify the molecule dictionary, and drop keys when the atoms are assigned to a new molecule
            while len(set(self.mol_dict.keys())) > self.nfrags:
                dist = distances()    
                dist.sort(key = lambda item: item[2]) # sort on distance of [i, j, dist]
                # take the second atom, append it to the list of atoms in the value of the first
                mol1, mol2 = dist[0][:2] # dist = [[i,j,dist], [i,j,dist], [i,j,dist]...]
                # adding all atoms in 'mol2' to the list of atoms of the first molecule
                for index, atom in enumerate(self.mol_dict[mol2]):
                    self.mol_dict[mol1].coords.append(atom)
                # remove the index of the molecule that has been added to the first list
                self.mol_dict.pop(mol2)   

        collect_atoms()
        # assign a number to each molecule
        for key, value in self.mol_dict.items():
            val = key
            for atom in value:
                atom.mol = val

    def check_db(self):
        """Checks fragments for a match in the database"""

        symbols = {} # not tied to an instance, only required when checking the database
        for frag, atoms in self.mol_dict.items():
            symbols[frag] = [atom.symbol for atom in atoms]
 
        def check_dict(molecules_dict, symbols_dict, db):
            """Checking molecule database for a match, and returning the required attributes- name, atoms in molecule, type of molecule, charge, multiplicity, and elements present in the molecule"""
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
                                "atoms": self.mol_dict[sym].coords, #not symbols, but the atom objects
                                "charge": charge,
                                "multiplicity": mult,
                                "elements": sort_elements(self.mol_dict[sym])
                            }
                        molecules_dict[sym] = data
            return molecules_dict

        self.fragments = {}
        for db in (Molecule.Anions, Molecule.Cations, Molecule.Neutrals):
            self.fragments = check_dict(self.fragments, symbols, db)

    def print_frags(self):
        for frag, data in self.fragments.items():
            print(f"{data['type'].capitalize()} found: {data['name']}")
    
    def separate(self):
        """Separates coordinates into specific fragments using the intermolecular distances. Note this function only works with intermolecular fragments and cannot split molecules on bonds."""
        self.split()
        self.check_db()
        self.print_frags() 
        self.fmo_meta() 

    #format for fmo run of complete system 
    def fmo_meta(self):
        """Creates strings for the INDAT and ICHARG blocks of GAMESS FMO calculations, bound to the
molecule instance as self.indat and self.charg"""
        
        info = {}
       
        #group together indat and charge for each fragment, and order according to the atom indices of indat 
        
        for frag, data in self.fragments.items():
            info[frag] = {"indat": f"0,{data['atoms'][0].index},-{data['atoms'][-1].index},",
                          "charg" : str(data['charge'])}
        # items need sorting
        # 0,1,7, ### sort on 2nd item ###
        # 0,8,28,
        # 0,29,35,
        # and not 
        # 0,1,7,
        # 0,29,35,
        # 0.8,28 
        
        sorted_info = sorted(info.items(), key = lambda val: int(val[1]['indat'].split(',')[1])) 
        # key: given an item (the val), sort by the value of the key, value pair (val[1]), then
        # select the indat, split on comma and return the second item, then convert to an integer
        # could also just sort on mol (or frag of info[frag]), from the assignments in self.split(), but these might not
        # always be in a numerical order- by using the index from self.coords, it is always ensured that the
        # correct order is shown, as these coords are also used in the input file        
        self.indat = []
        self.charg = []
        for val in sorted_info:
            self.indat.append(val[1]['indat'])
            self.charg.append(val[1]['charg'])
        self.indat.append('0')

def sort_elements(lst):
    els = []
    elements = set([atom.symbol for atom in lst])
    for i in elements:
        els.append((i, float(PT.get_atnum(i))))
    sorted_els = sorted(els, key = lambda val: val[1])
    return sorted_els
