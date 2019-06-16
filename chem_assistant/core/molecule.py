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

    Anions = {"bromide" : ["Br"]}
    Anions["chloride"] = ["Cl"]
    Anions["bf4"] = ['B', 'F', 'F', 'F', 'F']
    Anions["dca"] = ['N', 'C', 'N', 'C', 'N']
    Anions["pf6"] = ['F', 'P', 'F', 'F', 'F', 'F', 'F']
    Anions["mes"] = ['S', 'O', 'O', 'O', 'C', 'H', 'H', 'H']
    Anions["ntf2"] = ['F', 'F', 'F', 'F', 'F', 'N', 'S', 'S',
                        'O', 'O', 'O', 'O', 'C', 'C', 'F']
    Anions["bis-fsi"] = ['F', 'S', 'O', 'O', 'N', 'S', 'O', 'O', 'F']
    Anions["tos"] = ['C', 'C', 'C', 'C', 'H', 'H', 'H', 'H',
                        'H', 'H', 'H', 'S', 'O', 'O', 'O', 'C', 'C', 'C']
    Anions["dhp"] = ['H', 'H', 'P', 'O', 'O', 'O', 'O']
    Anions["acetate"] = ['C','H', 'H', 'H', 'C', 'O', 'O']
    Anions['saccharinate'] = ['C', 'C', 'C', 'C', 'C', 'C', 'H',
                     'H', 'H', 'H', 'C', 'O', 'N', 'S', 'O', 'O']

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
    Cations['lithium'] = ['Li']
    Cations['sodium'] = ['Na']
    Cations['potassium'] = ['K']

    Neutrals = {"nh3" : ['N', 'H', 'H', 'H']}
    Neutrals['water'] = ['H', 'H', 'O']
    Neutrals['dopamine-c=c-carbonyl'] = ['O','O','C','C','C','C','C','C','C','C','N','H','H','H','H','H']
    Neutrals['dopamine-c=c-hydroxyl'] = ['O','O','C','C','C','C','C','C','C','C','N','H','H','H','H','H','H','H']
    Neutrals['dopamine-c-c-hydroxyl'] = ['O','O','C','C','C','C','C','C','C','C','N','H','H','H','H','H','H','H','H','H']
    
    Radicals = {}


    def __init__(self, using = None, atoms = None, nfrags = None, user_created = False):
        if using is not None:
            self.coords = self.read_xyz(using)
            #list of Atom objects, more useful than list of coordinates
        self.bonds = []
        if atoms is not None and using is None:
            self.coords = atoms
        
        for index, atom in enumerate(self.coords):
            atom.index = index + 1 #so first atom has index 1
            atom.pos_in_list = index + 1

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
        return string[:-1]

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
        Note: Only works if all intramolecular bonds are shorter than all
        intermolecular bonds. A long bond in the middle of a molecule or a 
        short hydrogen bond will cause this to fail!
        """

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
            """Returns a list of lists of [i, j, min_dist(i, j)] for each combination of atoms in the system"""
            dist_list = []
            # here, decide if bonded!
            # for every combination of atoms in the list, find their distance
            for i, j in itertools.combinations(self.mol_dict.keys(), 2):
                minimum_dist = min_dist(self.mol_dict[i], self.mol_dict[j])
                dist_list.append([i, j, minimum_dist])
            return dist_list

        def collect_atoms():
            # modify the molecule dictionary, and drop keys when the atoms are assigned to a new molecule
            if self.nfrags is None:
                done = False
                while not done:
                    try:
                        self.nfrags = int(input('Number of fragments: '))
                        done = True
                    except ValueError:
                        print('Must give an integer number of fragments!')
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
            elif db == Molecule.Radicals:
                charge = 0
                mult = 2
                mol_type = 'radical'
            ## create more clauses for biradicals

            for name, atom_list in db.items():
                for sym, molecule in symbols_dict.items():
                    if sorted(molecule) == sorted(atom_list): #sorts in place, but no overwriting of variable
                        
                        data = {
                                "type": mol_type,
                                "name" : name,
                                "atoms": self.mol_dict[sym], #not symbols, but the atom objects
                                "charge": charge,
                                "multiplicity": mult,
                                "elements": sort_elements(self.mol_dict[sym]),
                                "frag_type": "frag"
                            }
                        molecules_dict[sym] = data
            return molecules_dict

        self.fragments = {}
        for db in (Molecule.Anions, Molecule.Cations, Molecule.Neutrals, Molecule.Radicals):
            self.fragments = check_dict(self.fragments, symbols, db)

        #sort order of atoms
        for data in self.fragments.values():
            data['atoms'] = sorted(data['atoms'], key = lambda atom: atom.index)
            # add number
            for i, atom in enumerate(data['atoms']):
                atom.number = i + 1


    def renumber_molecules(self):
        """
        Molecule numbers (Mol: _) are sometimes not in a numerical order. 
        This function takes the molecules and gives them a number from 1 to the number of fragments
        """ 

        current = set([atom.mol for atom in self.coords])
        convert_keys = {k: v for k, v in enumerate(self.fragments.keys(), 1)} # old: new
        frags = list(self.fragments.items())
        self.fragments.clear()
        for k, v in frags:
            for key, val in convert_keys.items():
                if k == val:
                    self.fragments[key] = v
                    for atom in self.fragments[key]['atoms']:
                        atom.mol = key

    def sort_fragments_by_index(self):
        """
        Sorts self.fragments according to the index of the atoms in each fragment.
        Sorts by the index of the first atom in each fragment
        """
        frags = [(k, v) for k,v in self.fragments.items()]
        frags = sorted(frags, key = lambda kv: kv[1]['atoms'][0].index)
        self.fragments = {k : v for k, v in frags}        
        
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

            manual_assignment = [i.strip() for i in manual_split]
            frag_indices = []

            start = 0
            end = 0
            # # pair up, but how if only one element?
            for frag in manual_assignment:
                if '[' in frag:
                    start = int(frag[1:])
                elif ']' in frag:
                    end = int(frag[:-1])
                if start != 0 and end != 0:
                    frag_indices.append((start, end))
                    start, end  = 0, 0 # start looking again
                elif start != 0 and end == 0:
                    # found first atom
                    pass
                else:
                    frag_indices.append(int(frag)) # single number

            return frag_indices 
            # [(start of frag, end of frag), single atom, (...), (...)]

        
        def update_mol_dictionary(frag_indices):
            for molecule, item in enumerate(frag_indices):
                num = molecule + 1
                if isinstance(item, tuple): # if molecule (1, 21) --> choline
                    start, end = item
                    for ind in range(start, end + 1):
                        self.coords[ind - 1].index = ind
                        self.coords[ind - 1].mol = num
                else: # single-atom 'molecule' i.e. Li, Na, K, Cl, Br
                    self.coords[item - 1].index = item
                    self.coords[item - 1].mol = num
        
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
2 fragments of water, type "[1, 3], [4, 6]", without quotes.
Give the atom numbers of the first and last atom in each fragment.

If you have ions of one element, say Lithium, included with the two water 
molecules, include the number without brackets: [1, 3], 4, [5, 7]
""")
        print()
        frag_indices = get_manual_assignments()
        update_mol_dictionary(frag_indices)
    
    def all_atoms_assigned(self):
        s = 0
        for frag in self.fragments.values():
            for atom in frag['atoms']:
                s += 1
        return s == len(self.coords)

    def assign_neighbours(self):
        """
        Checks each atom, either per fragment or in whole list, for bonded atoms by considering separation and van der waals radii
        """

        for i, atom_i in enumerate(self.coords):
            for j, atom_j in enumerate(self.coords):
                if i != j:
                    dist = atom_i.distance_to(atom_j)
                    vdw_dist = PT.get_vdw(atom_i) + PT.get_vdw(atom_j)
                    if dist < vdw_dist:
                        atom_i.connected_atoms.add(atom_j)
                        atom_j.connected_atoms.add(atom_i)
                        # also finds hydrogen bonds to atoms in other molecules

        
    def add_ionic_network(self):
        """Adds one item to self.fragments- removing all neutral species, along with the one-atom
        ions"""
        coord_list = [coord for coord in self.coords] 
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
            self.ionic = {
                "type": 'ionic',
                "name" : 'ionic',
                "atoms": coord_list, #not symbols, but the atom objects
                "charge": 0,
                "multiplicity": 1,
                "elements": sort_elements(coord_list),
                "frag_type": "ionic"
            }
   
    def separate_old(self):
        """
        Separates coordinates into specific fragments using the intermolecular distances. Note
        this function only works with intermolecular fragments and cannot split molecules on bonds.
        
        Also often couldn't determine fragments due to long intramolecular bonds. Now deprecated in
        favour of checking van der waals distances
        """
        self.split()
        self.check_db()
        self.renumber_molecules()
        self.print_frags() 
        all_assigned = False
        while not all_assigned:
            self.reassign_frags_manually()
            if self.all_atoms_assigned():
                all_assigned = True
        self.add_ionic_network()

    def separate(self):
        """
        Separates coordinates into specific fragments using the intermolecular distances along
        with van der waals radii. Note this function only works with intermolecular fragments 
        and cannot split molecules on bonds.
        """
        self.split_vdw()
        self.check_db()
        self.sort_fragments_by_index()
        self.renumber_molecules()
        # self.print_frags() 
        self.add_ionic_network()


    def distance_matrix(self):
        """
        Creates an N x N matrix of interatomic distances
        between every atom in the system. N = number of 
        atoms in system.
        """
        num_atoms = len(self.coords)
        matrix = np.zeros((num_atoms, num_atoms))

        for i, atom_i in enumerate(self.coords):
            for j, atom_j in enumerate(self.coords):
                    matrix[i, j] = atom_i.distance_to(atom_j)
        return matrix

    def split_vdw(self):

        mol_count = 0
        self.mol_dict = {}
        dists = self.distance_matrix()
        for i, atom_i in enumerate(self.coords):
            for j, atom_j in enumerate(self.coords):
                if i != j:
                    vdw_dist = PT.get_vdw(atom_i) + PT.get_vdw(atom_j)
                    if dists[i,j] < vdw_dist:                    
                        # bonded
                        # if atom_i already in dict, add j to same molecule
                        # or create new molecule
                        i_added = False
                        j_added = False

                        # does not prevent duplicates!
                        for k, v in self.mol_dict.items():
                            if atom_i in v and atom_j not in v:
                                v.append(atom_j)
                                j_added = True
                            if atom_j in v and atom_i not in v:
                                v.append(atom_i)
                                i_added = True

                        if not i_added and not j_added:
                            self.mol_dict[mol_count] = [atom_i, atom_j]
                            mol_count += 1
                    else:
                        # if non-bonded, check if already in dict
                        atom_i_in_dict = False
                        atom_j_in_dict = False

                        for k, v in self.mol_dict.items():
                            if atom_i in v:
                                atom_i_in_dict = True
                            if atom_j in v:
                                atom_j_in_dict = True
                        if not atom_i_in_dict:
                            self.mol_dict[mol_count] = [atom_i]
                            mol_count += 1
                        if not atom_j_in_dict:
                            self.mol_dict[mol_count] = [atom_j]
                            mol_count += 1

        # check if any index appears twice
        # if yes, delete it
        # if any other atoms are present, 
        # then add them to first molecule with index in
        # and remove from original place
        # should result in some empty lists, so delete them

        for k, v in self.mol_dict.items():
            indices_of_first_mol = set([atom.index for atom in v])
            for atom in v:
                # check all other 'molecules'
                for k2, v2 in self.mol_dict.items():
                    if k != k2:                        
                        for ind, atom2 in enumerate(v2):
                            if atom.index == atom2.index:
                                # if duplicate
                                del v2[ind]
                                # check all other atoms in second molecule
                                if atom2.index not in indices_of_first_mol:
                                    removed = v2.pop(ind)
                                    v.append(removed)
        
        # for some reason- still the odd atom left over...
        # check again for distances
        for index, mol1 in self.mol_dict.items():
            for index2, mol2 in self.mol_dict.items():
                if index != index2:
                    for atom1 in mol1:
                        for atom2 in mol2:
                            vdw_dist = PT.get_vdw(atom1) + PT.get_vdw(atom2)
                            if dists[atom1.index - 1, atom2.index - 1] < vdw_dist:
                                mol1.append(atom2)
                                mol2.remove(atom2)

        # remove empty lists, sort atoms by index
        self.mol_dict = {k:v for k,v in self.mol_dict.items() if len(v)>0}
        
        for mol in self.mol_dict.values():
            mol.sort(key = lambda atom: atom.index)

    def find_h_bonds(self):

        def find_partners(self):

            self.assign_neighbours()
            frag_list = [frag['atoms'] for frag in self.fragments.values()]

            h_bonders = ['O', 'F', 'H', 'N']
            H_BOND_DIST = 2.0
            
            # in one case a methyl hydrogen was 1.998Å away from O;
            # need to correct for that! Look at connecting atoms
            counted = []
            h_bonded = []
            for i, mol in enumerate(frag_list):
                for j, mol2 in enumerate(frag_list):
                    if i != j:
                        for atom1 in mol: # for every atom in first frag
                            for atom2 in mol2: # for every atom in second frag
                                if atom1.distance_to(atom2) < H_BOND_DIST:
                                    pairs = [(atom1.pos_in_list, atom1.mol), (atom2.pos_in_list, atom2.mol)]
                                    pairs.sort() #sort by pos_in_list first- doesn't matter how it sorts, only that it is applied to every atom in list


                                    # remove alkyl proton interactions

                                    if atom1.symbol == 'H' and 'C' in [atom.symbol for atom in atom1.connected_atoms]:
                                        pass
                                    if atom2.symbol == 'H' and 'C' in [atom.symbol for atom in atom2.connected_atoms]:
                                        pass

                                    if atom1.symbol and atom2.symbol in h_bonders and pairs not in counted:
                                        
                                        h_bonded.append([atom1, atom2, atom1.distance_to(atom2)])
                                        counted.append(pairs) 

                                        


                            # assigning twice- could be cut down
                            # could say for atom2 in mol2 if mol2 != mol- and
                            # change the iteration above?

            db = {}
            db['Cat'] = Molecule.Cations
            db['An']  = Molecule.Anions
            db['Neu'] = Molecule.Neutrals
            db['Rad'] = Molecule.Radicals

            # all combinations of cations/anions/radicals/neutrals
            combinations = list(itertools.combinations_with_replacement(db.keys(), 2))              

            self.hbond_data_to_export = []            

            for bond in h_bonded:
                one, two, dist = bond
                one_name = self.fragments[one.mol]['name']
                two_name = self.fragments[two.mol]['name']
                # to form a h-bond, must contain two of the h-bonding atoms in that molecule (O-H bonds, N-H bonds)
 
                for comb in combinations:
                    group1, group2 = comb

                    if any(name in db[group1] for name in (one_name, two_name)) and any(name in db[group2] for name in (one_name, two_name)):
                        print(f'{group1}-{group2}')
                        print(f"({one.symbol}, {self.fragments[one.mol]['name']} (mol {one.mol}))--- {dist:.3f} Å ---({two.symbol}, {self.fragments[two.mol]['name']} (mol {two.mol}))")
                        self.hbond_data_to_export.append([self.fragments[one.mol]['name'], 
                                           one.symbol, 
                                           self.fragments[two.mol]['name'], 
                                           two.symbol, 
                                           dist])

            

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
                        
        def bond_components(self, tup):
            one, two, length = tup
            atom = self.coords[one - 1]
            atom2 = self.coords[two - 1]    

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

        find_partners(self)
        return self.hbond_data_to_export
