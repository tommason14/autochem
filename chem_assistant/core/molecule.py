#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .periodic_table import PeriodicTable as PT
from .atom import Atom
from .bond import Bond
from .utils import sort_elements

import re
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
    Anions["dhp"] = ['H', 'H', 'P', 'O', 'O', 'O', 'O']
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


    def __init__(self, using = None, atoms = None, nfrags = None, user_created = False):
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
                                "elements": sort_elements(self.mol_dict[sym]),
                                "frag_type": "frag"
                            }
                        molecules_dict[sym] = data
            return molecules_dict

        self.fragments = {}
        for db in (Molecule.Anions, Molecule.Cations, Molecule.Neutrals):
            self.fragments = check_dict(self.fragments, symbols, db)

        #sort order of atoms
        for data in self.fragments.values():
            data['atoms'] = sorted(data['atoms'], key = lambda atom: atom.index)
            # add number
            for i, atom in enumerate(data['atoms']):
                atom.number = i + 1


    def renumber_molecules(self):
        """Molecule numbers (Mol: _) are sometimes not in a numerical order. This function takes the
        molecules and gives them a number from 1 to the number of fragments""" 
        current = set([atom.mol for atom in self.coords])
        converter = {k: v for k, v in enumerate(current, 1)}

        convert_keys = {k: v for k, v in enumerate(self.fragments.keys(), 1)} # old: new
        frags = list(self.fragments.items())
        self.fragments.clear()
        for k, v in frags:
            for key, val in convert_keys.items():
                if k == val:
                    self.fragments[key] = v
       
        for atom in self.coords:
            for k, v in converter.items():
                if atom.mol == v:
                    atom.mol = k
        
    def print_frags(self):
        print()
        for frag, data in self.fragments.items():
            print(f"{data['type'].capitalize()} found: {data['name']}")
        print()

    def reassign_frags_manually(self):
        """Called if fragments are not assigned correctly- user then input the fragments manually"""

        def get_manual_assignments():
            manual = input("Fragments: ")
            manual_split = manual.split(',')
            manual_assignment = [int(i.strip()) for i in manual_split]
            frag_indices = []
            for i in range(0, len(manual_assignment) - 1, 2):
                frag_indices.append((manual_assignment[i], manual_assignment[i + 1]))
            return frag_indices
        
        def update_mol_dictionary(frag_indices):
            for molecule, pair in enumerate(frag_indices):
                num = molecule + 1
                start, end = pair
                for ind in range(start, end + 1):
                    self.coords[ind - 1].index = ind
                    self.coords[ind - 1].mol = num
            # reassign mol dict, check db
            self.mol_dict.clear()
            mols = set([atom.mol for atom in self.coords])
            for mol in mols:
                self.mol_dict[mol] = Molecule(atoms = [atom for atom in self.coords if atom.mol == mol])
            self.check_db()
            self.print_frags()

        print()
        print(f"Should have found {self.nfrags} fragments, but found {len(self.fragments)}...")
        print()
        print("""\
Fragments are too close together- there are bonds within fragments
that are longer than the shortest distance between fragments.

Type in the fragments manually. For example, if you have 
2 fragments of water, type "1, 3, 4, 6", without quotes.
Give the atom numbers of the first and last atom in each fragment""")
        print()
        frag_indices = get_manual_assignments()
        update_mol_dictionary(frag_indices)
    
    def all_atoms_assigned(self):
        s = 0
        for frag in self.fragments.values():
            for atom in frag['atoms']:
                s += 1
        return s == len(self.coords)

        
    def add_ionic_network(self):
        """Adds one item to self.fragments- removing all neutral species, along with the one-atom
ions"""
        coord_list = [coord for coord in self.coords] # copy of the list
        for k, frag in self.fragments.items():
            # remove neutrals
            if frag['charge'] is 0:
                for coord in frag['atoms']:
                    # remove coord from coord_list
                    for i, atom in enumerate(coord_list):
                        if coord.index == atom.index:
                            del coord_list[i]
            # remove Li, Na, Cl, Br etc...
            elif len(frag['atoms']) == 1:
                for coord in frag['atoms']:
                    # remove coord from coord_list
                    for i, atom in enumerate(coord_list):
                        if coord.index == atom.index:
                            del coord_list[i]
        if len(coord_list) == len(self.coords):
            pass
        else:
            self.fragments[len(self.fragments) + 1] = {
                "type": 'ionic',
                "name" : 'ionic',
                "atoms": coord_list, #not symbols, but the atom objects
                "charge": 0,
                "multiplicity": 1,
                "elements": sort_elements(coord_list),
                "frag_type": "ionic"
            }
   
    def separate(self):
        """Separates coordinates into specific fragments using the intermolecular distances. Note this function only works with intermolecular fragments and cannot split molecules on bonds."""
        self.split()
        self.check_db()
        self.renumber_molecules()
        self.print_frags() 
        if not self.all_atoms_assigned():
            self.reassign_frags_manually()
        self.add_ionic_network()
        
    def find_h_bonds(self):
  
        frag_list = [frag['atoms'] for frag in self.fragments.values()]

        h_bonders = ['O', 'F', 'H', 'N']
        H_BOND_DIST = 2.0 
        
        # in one case a methyl hydrogen was 1.998Å away from O;
        # need to correct for that! Look at connecting atoms

        for i, mol in enumerate(frag_list):
            for j, mol2 in enumerate(frag_list):
                if i != j:
                    for atom1 in mol: # for every atom in first frag
                        for atom2 in mol2: # for every atom in second frag
                            if atom1.distance_to(atom2) < H_BOND_DIST:
                                if atom1.symbol and atom2.symbol in h_bonders and atom1.symbol is not atom2.symbol:
                                    atom1.h_bonded_to.append(atom2) 
                        # assigning twice- could be cut down
                        # could say for atom2 in mol2 if mol2 != mol- and change the iteration above?


        def remove_duplicate_bonds(self):
            h_bonded = [] #entire system
            for atom in self.coords:
                if len(atom.h_bonded_to) > 0:
                    per_atom = [] # each atom- may have more than one h-bonding partner
                    for atom2 in atom.h_bonded_to:
                        per_atom.append((atom.index, atom2.index, atom.distance_to(atom2)))    
                    h_bonded.append(per_atom)

            # now should flatten list- easier to remove duplicates that way
            flattened = []
            for lst in h_bonded:
                for l in lst:
                    flattened.append(l)

            # remove duplicate h-bonds
            # make copy first:
            h_bonded_flat = flattened
            for i, bond in enumerate(h_bonded_flat):
                for j, bonds in enumerate(h_bonded_flat):
                    if i != j:
                        if sorted(bond[:2]) == sorted(bonds[:2]): # if same atoms are involved
                            del h_bonded_flat[j] # remove second instance
            return h_bonded_flat

        hbonds = remove_duplicate_bonds(self)
                        
        def bond_components(self, tup):
            one, two, length = tup
            atom = self.coords[one - 1]
            atom2 = self.coords[two - 1]    

            # print(f"({atom.symbol}, mol: {atom.mol}, atom: {atom.number})--- {length:.3f}Å ---({atom2.symbol}, mol: {atom2.mol}, atom: {atom2.number})")

            mol1 = self.fragments[atom.mol]['name']
            mol2 = self.fragments[atom2.mol]['name']
        
            print(f"({atom.symbol}, {mol1} (mol {atom.mol}))--- {length:.3f}Å ---({atom2.symbol}, {mol2} (mol {atom2.mol}))")

            # Choline-Acetate rather than Acetate-Choline
            if mol1.lower() in Molecule.Anions:
                connection = f"{mol2}-{mol1}"
            else:
                connection = f"{mol1}-{mol2}"
            length = round(length, 3)
            return [connection, length]

        bonds_in_mol = []
        print('Hydrogen bonds:')
        for bond in hbonds:
            data = bond_components(self, bond)
            bonds_in_mol.append(data) # list, so can add file and path later
        print() # nicer formatting from command line
        return bonds_in_mol
