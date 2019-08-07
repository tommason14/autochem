from .periodic_table import PeriodicTable as PT
from .atom import Atom
from .utils import sort_elements

import re
import numpy as np
import math
import itertools
import sys

__all__ = ['Molecule']

class Molecule:
    """
    Implementing the concept of a molecular system. Crucial is the
    separation of molecules in a given system, with output for Fragment
    Molecular Orbital calculations
    """

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
    Anions['triflate'] = ['C', 'F', 'F', 'F', 'S', 'O',  'O', 'O']

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
    Neutrals['benzene'] = ['C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H'] 
    Neutrals['acetone'] = ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'O']
    Neutrals['dmso'] = ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'S', 'O']
    
    Radicals = {}


    def __init__(self, using = None, atoms = None, nfrags = None, user_created = False):
        if using is not None:
            self.coords = self.read_xyz(using)
        if atoms is not None and using is None:
            if len(atoms) == 0:
                sys.exit('Error: atoms argument passed into Molecule is empty')
            if not isinstance(atoms[0], Atom):
                self.coords = []
                for atom in atoms:
                    sym, *coords = atom
                    a = Atom(symbol = sym, coords = coords)
                    self.coords.append(a)
            else: 
                self.coords = atoms
    
        for index, atom in enumerate(self.coords):
            atom.index = index + 1

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
        """
        Repr for system after fragmentation
        """
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
        """
        Returns the molecular format in a variety of formats using keyword 
        arguments:
        * *as_dict* -- returns a python dictionary
        * *as_latex* -- returns a latex expression using \textsubscript{} notation (equivalent to $_{value}$, but the math expression renders in a different font)
        * *as_html* -- returns a html expression using <sub> tags

        If no keyword arguments are passed, a regular string is returned with 
        no formatting i.e. C8H18
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
        """
        Returns molecular mass in g mol⁻¹
        """
        formula = self.formula(as_dict=True)
        mass = 0
        for element, number in formula.items():
            mass += PT.get_mass(element) * number
        return f"{mass:.2f} g mol⁻¹"

    def overall_charge_and_mult(self):
        """
        Checks system for overall charge and multiplicity
        """
        mol = Molecule(atoms = [i for i in self.coords])
        mol.separate()
        charge = Molecule.get_charge(mol.fragments)
        multiplicity = Molecule.get_multiplicity(mol.fragments)
        print(charge, multiplicity)
        return charge, multiplicity
     
    def read_xyz(self, using):
        """
        Reads coordinates of an xyz file and return a list of |Atom| objects,
        one for each atom
        """
        coords = []
        with open(using, "r") as f:
            for coord in f.readlines()[2:]:
                line = coord.split()
                for val in PT.ptable.values():
                    if line[0] == val[0]:
                        coords.append(Atom(line[0], coords = tuple(float(i) for i in line[1:4])))
        return coords

    def write_xyz(self, atoms, filename = None):
        """
        Writes an xyz file using a list of |Atom| instances
        """
        if filename is None:
            raise ValueError('write_xyz: Must give a path to the output file')
        else:
            with open(filename, "w") as f:
                f.write(str(len(atoms)) + '\n\n')
                for atom in atoms:
                    f.write(f"{atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f} \n")


    def split_old(self):
        """
        Assigns each atom in ``self.coords`` to a different fragment- very 
        useful for FMO calculations.
        Note: Only works if all intramolecular bonds are shorter than all
        intermolecular bonds. A long bond in the middle of a molecule or a 
        short hydrogen bond will cause this to fail!
        """

        self.mol_dict = {}
        for atom in self.coords:
            self.mol_dict[atom.index] = Molecule(atoms = [atom]) 
            # 1-atom molecules

        def min_dist(mol_i, mol_j):
            """
            Finds the minimum distance between atoms in two different molecules
            """
            distances = []
            for atom_i in mol_i:
                for atom_j in mol_j:
                    distances.append(atom_i.distance_to(atom_j))
            return min(distances)

        def distances():
            """
            Returns a list of lists of [i, j, min_dist(i, j)] for each
            combination of atoms in the system
            """
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
                dist.sort(key = lambda item: item[2]) # sort on distance                                 
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
        """
        Checks fragments for a match in the database
        """

        symbols = {} 
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
                    if sorted(molecule) == sorted(atom_list): 
                        data = {
                                "type": mol_type,
                                "name" : name,
                                "atoms": self.mol_dict[sym],
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
        This function takes the molecules and gives them a number from 1 to the 
        number of fragments
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
        """
        Called if fragments are not assigned correctly- user then inputs 
        the fragments manually
        """

        def get_manual_assignments():
            manual = input("Fragments: ")
            manual_split = manual.split(',')
            manual_assignment = [i.strip() for i in manual_split]
            frag_indices = []

            start = 0
            end = 0
            # pair up, but how if only one element?
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

    def assign_neighbours(self):
        """
        Checks each atom, either per fragment or in whole list, for bonded 
        atoms by considering separation and van der waals radii
        """

        for i, atom_i in enumerate(self.coords):
            for j, atom_j in enumerate(self.coords):
                if i != j:
                    dist = atom_i.distance_to(atom_j)
                    vdw_dist = PT.get_vdw(atom_i) + PT.get_vdw(atom_j)
                    if dist < vdw_dist:
                        if atom_j not in atom_i.connected_atoms:
                            atom_i.connected_atoms.append(atom_j)
                        if atom_i not in atom_j.connected_atoms:
                            atom_j.connected_atoms.append(atom_i)


    def add_ionic_network(self):
        """
        Adds one item to self.fragments- removing all neutral species, 
        along with the one-atom ions
        """

        def ionic_mol_properties(coords):
            """
            Returns charge and multiplicity of ionic network
            """
            ionic_mol = Molecule(atoms = coords)
            # causing errors here I think, when separating and using for input
            ionic_mol.separate()
            ionic_frags = ionic_mol.fragments
            charge = sum(frag['charge'] for frag in ionic_frags.values())
            multiplicity = Molecule.get_multiplicity(ionic_frags)

            return charge, multiplicity

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
        if len(coord_list) != len(self.coords) and len(coord_list) != 0:
            # split and add charges and multiplicities up
            charge, multiplicity = ionic_mol_properties(coord_list)
            self.ionic = {
                "type": 'ionic',
                "name" : 'ionic',
                "atoms": coord_list,
                "charge": charge,
                "multiplicity": multiplicity,
                "elements": sort_elements(coord_list),
                "frag_type": "ionic"
            }

    def all_atoms_assigned(self):
        """
        Checks self.coords to see if all atoms have a molecule assigned
        """
        for atom in self.coords:
            if atom.mol is None:
                return False
        return True
   
    def separate(self):
        """
        Separates coordinates into specific fragments using the intermolecular 
        distances along with van der waals radii. Note this function only works 
        with intermolecular fragments and cannot split molecules on bonds.
        """

        self.split()
        self.check_db()
        self.renumber_molecules()
        self.sort_fragments_by_index()
        if not self.all_atoms_assigned():
            print(f"{len(self.fragments)} fragments found, some atoms unaccounted for...")
            all_assigned = False
            while not all_assigned:
                self.reassign_frags_manually()
                if self.all_atoms_assigned():
                    all_assigned = True
        # self.add_ionic_network()


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

    def split(self):
        """
        Split a system into fragments using van der waals radii. Modifies
        attributes of atoms in self.coords directly in loop, instead of creating
        a dictionary and appending to the dictionary as we go. This
        significantly speeds up the fragmentation.
        """
        mol_count = 0
        dists = self.distance_matrix()
        for i, atom_i in enumerate(self.coords):
            connected = False
            for j, atom_j in enumerate(self.coords):
                if i != j:
                    vdw_dist = PT.get_vdw(atom_i) + PT.get_vdw(atom_j)
                    if dists[i,j] < vdw_dist:
                        # connected
                        connected = True
                        if atom_i not in atom_j.connected_atoms:
                            atom_j.connected_atoms.append(atom_i)
                        if atom_j not in atom_i.connected_atoms:
                            atom_i.connected_atoms.append(atom_j)
                        if atom_i.mol is None and atom_j.mol is None:
                            atom_i.mol = mol_count
                            atom_j.mol = mol_count
                            mol_count += 1
                        elif atom_i.mol is None and atom_j.mol is not None:
                            atom_i.mol = atom_j.mol
                        elif atom_j.mol is None and atom_i.mol is not None:
                            atom_j.mol = atom_i.mol
                        # if different assignments, remove original assignment
                        # combine the two fragments together, as they are connected
                        elif atom_i.mol is not None and atom_j.mol is not None:
                            if atom_i.mol != atom_j.mol:
                                orig = atom_j.mol
                                for atom in self.coords:
                                    if atom.mol is orig:
                                        atom.mol = atom_i.mol
            if not connected:
                atom_i.mol = mol_count
                mol_count += 1

        nums = set([atom.mol for atom in self.coords])
        self.mol_dict = {val: [atom for atom in self.coords if atom.mol == val] for val in nums}
        
        for mol in self.mol_dict.values():
            mol.sort(key = lambda atom: atom.index)
        
    @classmethod
    def create_complex_properties_dict(self):
        """
        Checks system for overall charge and multiplicity, and returns
        properties in a dictionary
        """
        def check_charge_and_mult(atoms):
            """
            Must check charge of overall complex, may not be zero.
            This function achieves this by creating a copy of the molecule, 
            fragmenting the copy and then summing the charges of each fragment.
            The reason being that the user may not want to have each fragment
            written to a file- which would happen if fragmenting self??? Test it
            Modify the self.complex dictionary
            """
            mol = Molecule(atoms = atoms)
            mol.separate()
            charge = Molecule.get_charge(mol.fragments)
            multiplicity = Molecule.get_multiplicity(mol.fragments)
            print(charge, multiplicity)
            return charge, multiplicity

        charge, mult = check_charge_and_mult(atoms)
        # self.complex used in input files
        # assuming a neutral closed shell system
        return {
            "type": "complex",
            "name": "complex",
            "atoms": atoms,
            "charge": charge,
            "mult": mult,
            "elements": sort_elements(atoms)
        }


    def find_h_bonds(self):
        """
        Gives hydrogen bonding data back to the user. Checks for suitable
        connected atoms and bond lengths of less than 2 Å, and bond angles of
        45° either side of linear.
        """

        def find_bonds(self):


            def valid_atoms(atom1, atom2):
                """
                Checks if a bond is formed between a pair of atoms of the
                correct type. This is achieved by checking the atoms
                that the bonds are connected to, not involved in the bond.
                For example, if an alkyl chain is close to an anion, this
                function ignores any of those interactions.
                """
                # remove alkyl proton interactions
                # need the atom.symbol == line, 
                # otherwise c-o --- h-o hydrogen
                # bonds will fail (as c is not in 
                # h-bonders)
                h_bonders = ['O', 'F', 'H', 'N']

                for atom in (atom1, atom2):
                    if atom.symbol not in h_bonders:
                        return False
                    if atom.symbol == 'H':
                        for a in atom.connected_atoms:
                            if a.symbol not in h_bonders:
                                return False
                return True

            def within_hbond_distance(atom1, atom2):
                """
                Checks that atoms are within hydrogen-bonding distances, set to
                2 Å.
                """
                H_BOND_DIST = 2.0

                return atom1.distance_to(atom2) < H_BOND_DIST

            def bond_angle(atom1, atom2):
                """
                Returns the angle between two atoms.

                 A
                 \ 
                   B --- C

                Checks angle ∠ABC, when A is atom1, B is atom2 and C is the
                connecting atom to atom2. There will probably be only one
                connecting atom, but if there is more, only considering the
                first in the list.
                """
                connected_to_atom2 = atom2.connected_atoms[0]
                return atom2.angle_between(atom1.coords,connected_to_atom2.coords)
                

            def within_angle_tolerance(atom1, atom2):
                """
                Checks that a hydrogen bond is formed at a suitable angle, 45°
                either side of linear.
                """
                return 225 > bond_angle(atom1, atom2) > 145


            def valid_bond(atom1, atom2):
                """
                Checks that two atoms forms a valid hydrogen bond.
                """
                if within_hbond_distance(atom1, atom2) and \
                within_angle_tolerance(atom1, atom2) and \
                valid_atoms(atom1, atom2):
                    return True
                return False

            self.assign_neighbours()
            frag_list = [frag['atoms'] for frag in self.fragments.values()]

            counted = []
            h_bonded = []
            for i, mol in enumerate(frag_list):
                for j, mol2 in enumerate(frag_list):
                    if i != j:
                        for atom1 in mol: 
                            for atom2 in mol2: 
                                if valid_bond(atom1, atom2):
                                    pair = [atom1.index, atom2.index]
                                    pair.sort() 
                                    if pair not in counted:
                                        dist = atom1.distance_to(atom2)
                                        angle = bond_angle(atom1, atom2)
                                        h_bonded.append([atom1, atom2, dist, angle])
                                        counted.append(pair) 

            return h_bonded

        def find_molecule_type(molecule):
            """
            Finds the type of molecule in terms of charge and multiplicity.
            Returns 'Cat', 'An', 'Neu' or 'Rad'. Note that radicals are doublets
            of zero charge.
            """
            db = {}
            db['Cat'] = Molecule.Cations
            db['An']  = Molecule.Anions
            db['Neu'] = Molecule.Neutrals
            db['Rad'] = Molecule.Radicals

            for moltype, molecules in db.items():
                if molecule in molecules:
                    return moltype


        def hydrogen_bond_data(self, h_bonds):
            """
            Prints data to the screen, and returns a list of hydrogen bond
            attributes of the form:
            molecule1, atom1, molecule2, atom2, distance.
            """
            hbond_data = []            
            

            for bond in h_bonds:
                one, two, dist, angle = bond
                mol_one_name = self.fragments[one.mol]['name']
                mol_two_name = self.fragments[two.mol]['name']
            
                group1 = find_molecule_type(mol_one_name)
                group2 = find_molecule_type(mol_two_name)

                hbond_data.append([
                    mol_one_name, 
                    one.symbol, 
                    mol_two_name, 
                    two.symbol, 
                    dist,
                    angle
                ])

            return hbond_data

        hbonds = find_bonds(self)
        self.hbond_data = hydrogen_bond_data(self, hbonds)

        return self.hbond_data

    def frag_name(self, atom):
        """
        Accepts an atom from the molecule, returning the name of the fragment containing that atom.
        """

        if not hasattr(self, 'fragments'):
            self.separate()
        for frag in self.fragments.values():
            for a in frag['atoms']:
                if a.index == atom.index:
                    return frag['name']

    @classmethod
    def get_charge(cls, fragment_dict):
        return sum(frag['charge'] for frag in fragment_dict.values())

    @classmethod
    def get_multiplicity(cls, fragment_dict):
        return 2 if any(frag['type'] == 'radical' for frag in fragment_dict.values()) else 1
        # extend multiplicity for biradicals etc...

