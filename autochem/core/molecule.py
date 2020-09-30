from .periodic_table import PeriodicTable as PT
from .atom import Atom
from .utils import sort_elements

import re
import os
import numpy as np
import math
import itertools
import sys

__all__ = ['Molecule']


class Molecule:
    """
    Implementing the concept of a molecular system. The separation of molecules
    in a given system is also accounted for, important for Fragment
    Molecular Orbital calculations.

    Class Attributes
    ----------
    Anions: dict
        format of {'name': [atomic symbols]} for anions
    Cations: dict
        format of {'name': [atomic symbols]} for cations
    Neutrals: dict
        format of {'name': [atomic symbols]} for neutral molecules
    Radicals: dict
        format of {'name': [atomic symbols]} for doublet radicals
    Dications: dict
        format of {'name': [atomic symbols]} for ions of +2 charge
    Anion_radicals: dict
        format of {'name': [atomic symbols]} for negatively charged radicals
    Cation_radicals: dict
        format of {'name': [atomic symbols]} for positively charged radicals
    Dication_radicals: dict
        format of {'name': [atomic symbols]} for doubly charged radicals

    Instance Attributes
    -------------------
    xyz: string
        name of xyz file used to create the molecule
    coords: list 
        list of `Atom` objects for every atom in the molecule
    fragments: dict
        format of {number: subdict} created when `self.separate()` is called.
        The subdict contains the keys: type (string), name (string),
        atoms (list of `Atom` instances), charge (int), mult (int), 
        elements (list of atomic symbols)

    """

    Anions = {"bromide": ["Br"]}
    Anions["chloride"] = ["Cl"]
    Anions["bf4"] = ['B', 'F', 'F', 'F', 'F']
    Anions["dca"] = ['N', 'C', 'N', 'C', 'N']
    Anions["pf6"] = ['F', 'P', 'F', 'F', 'F', 'F', 'F']
    Anions["mes"] = ['S', 'O', 'O', 'O', 'C', 'H', 'H', 'H']
    Anions["ntf2"] = [
        'F', 'F', 'F', 'F', 'F', 'N', 'S', 'S', 'O', 'O', 'O', 'O', 'C', 'C',
        'F'
    ]
    Anions["bis-fsi"] = ['F', 'S', 'O', 'O', 'N', 'S', 'O', 'O', 'F']
    Anions["tos"] = [
        'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'S', 'O', 'O',
        'O', 'C', 'C', 'C'
    ]
    Anions["dhp"] = ['H', 'H', 'P', 'O', 'O', 'O', 'O']
    Anions["acetate"] = ['C', 'H', 'H', 'H', 'C', 'O', 'O']
    Anions['saccharinate'] = [
        'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'C', 'O', 'N', 'S',
        'O', 'O'
    ]
    Anions['triflate'] = ['C', 'F', 'F', 'F', 'S', 'O', 'O', 'O']

    Cations = {
        "c1mim": [
            'C', 'N', 'C', 'N', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H',
            'H', 'H', 'H'
        ]
    }
    Cations["c1mpyr"] = [
        'C', 'C', 'C', 'N', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H', 'H', 'H'
    ]

    Cations["c2mim"] = [
        'C', 'N', 'C', 'N', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H'
    ]
    Cations["c2mpyr"] = [
        'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'
    ]

    Cations["c3mim"] = [
        'N', 'C', 'N', 'C', 'C', 'C', 'C', 'H', 'C', 'C', 'H', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'
    ]
    Cations["c3mpyr"] = [
        'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'
    ]

    Cations["c4mim"] = [
        'C', 'N', 'C', 'C', 'N', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'
    ]
    Cations["c4mpyr"] = [
        'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
        'H', 'H'
    ]
    Cations['choline'] = [
        'N', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C',
        'H', 'H', 'C', 'H', 'H', 'O', 'H'
    ]
    Cations['lithium'] = ['Li']
    Cations['sodium'] = ['Na']
    Cations['potassium'] = ['K']
    Cations['styrene-trimethylammonium'] = [
        'H', 'C', 'H', 'C', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'C', 'H', 'C',
        'H', 'C', 'H', 'H', 'N', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C',
        'H', 'H', 'H'
    ]
    Cations['maotmac'] = [
        'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'C', 'O', 'O', 'C', 'H', 'H',
        'C', 'H', 'H', 'N', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'H',
        'H', 'H'
    ]

    Neutrals = {"nh3": ['N', 'H', 'H', 'H']}
    Neutrals['hf'] = ['H', 'F']
    Neutrals['methane'] = ['C', 'H', 'H', 'H', 'H']
    Neutrals['ethane'] = ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']
    Neutrals['water'] = ['H', 'H', 'O']
    Neutrals['acetic_acid'] = ['C', 'H', 'H', 'H', 'C', 'O', 'O', 'H']
    Neutrals['benzene'] = [
        'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H'
    ]
    Neutrals['acetone'] = ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'O']
    Neutrals['dmso'] = ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'S', 'O']
    Neutrals['amps'] = [
        'H', 'C', 'H', 'C', 'H', 'C', 'O', 'N', 'H', 'C', 'C', 'H', 'H', 'H',
        'C', 'H', 'H', 'H', 'C', 'H', 'H', 'S', 'O', 'H', 'O', 'O'
    ]
    Neutrals['sty-sulfonate-hydrogenated'] = [
        'H', 'C', 'C', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'C', 'H', 'C',
        'H', 'S', 'O', 'H', 'O', 'O'
    ]

    # Dopamine "monomers"
    Neutrals['dhica'] = [
        'O', 'C', 'C', 'C', 'N', 'H', 'C', 'C', 'C', 'C', 'C', 'O', 'H', 'H',
        'H', 'H', 'C', 'H', 'O', 'O', 'H'
    ]
    Neutrals['indole-5,6-dione'] = [
        'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'H', 'H', 'H',
        'H', 'H'
    ]
    Neutrals['dhi'] = [
        'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'H', 'H', 'H',
        'H', 'H', 'H', 'H'
    ]
    Neutrals['reduced-dhi'] = [
        'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H', 'H'
    ]
    Neutrals['dop-cov-dimer'] = [
        'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'N',
        'C', 'C', 'C', 'C', 'O', 'H', 'O', 'O', 'H', 'O', 'H', 'H', 'H', 'H',
        'H', 'H', 'H', 'H'
    ]
    Neutrals['uncyclised-dopamine'] = [
        'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'N',
        'H', 'H', 'H', 'H', 'O', 'H', 'O', 'H'
    ]

    Radicals = {
        'amps-dimer': [
            'H', 'O', 'S', 'O', 'O', 'C', 'H', 'H', 'C', 'C', 'H', 'H', 'H',
            'C', 'H', 'H', 'H', 'N', 'H', 'C', 'O', 'C', 'H', 'H', 'C', 'H',
            'H', 'C', 'H', 'H', 'C', 'H', 'C', 'O', 'N', 'H', 'C', 'C', 'H',
            'H', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'H', 'S', 'O', 'O', 'O',
            'H'
        ]
    }
    Radicals['divinyl-benzene'] = [
        'C', 'C', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C',
        'H', 'C', 'C', 'H', 'C', 'H', 'H'
    ]
    Radicals['ethylene-glycol-dimethacrylate'] = [
        'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'O', 'O', 'C', 'H',
        'H', 'C', 'H', 'H', 'O', 'C', 'O', 'C', 'C', 'H', 'H', 'C', 'H', 'H',
        'H'
    ]
    Radicals['glycerol-dimethacrylate'] = [
        'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'C', 'O', 'O', 'C', 'H',
        'H', 'C', 'H', 'O', 'H', 'C', 'H', 'H', 'O', 'C', 'O', 'C', 'C', 'H',
        'H', 'C', 'H', 'H', 'H'
    ]
    Radicals['sty-sulf-dimer-hydrogenated'] = [
        'H', 'O', 'S', 'O', 'O', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H',
        'C', 'C', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'C', 'H', 'C', 'C', 'H',
        'C', 'H', 'C', 'H', 'C', 'H', 'C', 'S', 'O', 'H', 'O', 'O', 'H'
    ]
    Dications = dict()
    Anion_radicals = dict()
    Cation_radicals = dict()

    Dication_radicals = {
        'maotmac-dimer-hydrogenated-rad': [
            'C', 'N', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H',
            'C', 'H', 'H', 'C', 'H', 'H', 'O', 'C', 'O', 'C', 'C', 'H', 'H',
            'H', 'C', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H',
            'C', 'O', 'O', 'C', 'H', 'H', 'C', 'H', 'H', 'N', 'C', 'H', 'H',
            'H', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H'
        ]
    }

    molecules = {
        **Cations,
        **Anions,
        **Neutrals,
        **Radicals,
        **Dications,
        **Anion_radicals,
        **Cation_radicals,
        **Dication_radicals
    }

    def __init__(self,
                 using=None,
                 atoms=None,
                 group=None,
                 bonds_to_split=None):
        self.check_user_additions()
        if using is not None:
            self.xyz = using
            self.coords = self.read_xyz(self.xyz)
        if atoms is not None and using is None:
            if len(atoms) == 0:
                sys.exit('Error: atoms argument passed into Molecule is empty')
            if not isinstance(atoms[0], Atom):
                self.coords = []
                for atom in atoms:
                    sym, *coords = atom
                    a = Atom(symbol=sym, coords=coords)
                    self.coords.append(a)
            else:
                self.coords = atoms

        self.frags_grouped_if_desired = False
        if group is not None:
            self.group_together = group

        self.split_on_bonds = False
        if bonds_to_split is not None:
            self.bonds_to_split = bonds_to_split
            self.split_on_bonds = True

        for index, atom in enumerate(self.coords):
            atom.index = index + 1

        if hasattr(self, 'coords'):
            # self.complex used in input files
            # assuming a neutral closed shell system as the default
            self.complex = {
                "type": "complex",
                "name": "complex",
                "atoms": self.coords,
                "charge": 0,
                "mult": 1,
                "elements": sort_elements(self.coords)
            }
        self.calc_overall_charge_and_mult()

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

    def translate(self, vector, frag=None):
        """
        Apply the vector to every atom in the system.
        Note that if fragmented, can specify which fragment to translate,
        by specifying a key of self.fragments
        """
        for atom in self.coords:
            atom.translate(vector)

        if frag is not None:
            if not hasattr(self, 'fragments'):
                raise AttributeError('Must run self.separate() first')
            for atom in self.fragments[frag]:
                atom.translate(vector)

    def formula(self, as_dict=False, as_latex=False, as_html=False):
        """
        Returns the molecular format in a variety of formats using keyword 
        arguments:
        * *as_dict* -- returns a python dictionary
        * *as_latex* -- returns a latex expression using \textsubscript{} notation (equivalent to $_{value}$, but the math expression renders in a different font)
        * *as_html* -- returns a html expression using <sub> tags

        If no keyword arguments are passed, a regular string is returned with 
        no formatting i.e. C8H18
        """

        formula = {}
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

    @property
    def mass(self):
        """
        Returns molecular mass in g mol⁻¹
        """
        formula = self.formula(as_dict=True)
        mass = 0
        for element, number in formula.items():
            element = Atom(element)  # PT.x requires Atom objects
            mass += PT.get_mass(element) * number
        return f"{mass:.2f} g mol⁻¹"

    def calc_overall_charge_and_mult(self):
        """
        Checks system for overall charge and multiplicity
        """
        if not hasattr(self, 'fragments'):
            self.separate()
        if self.split_on_bonds:
            self.fragment_on_bonds()
        self.overall_charge = Molecule.get_charge(self.fragments)
        self.overall_mult = Molecule.get_multiplicity(self.fragments)

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
                        coords.append(
                            Atom(line[0],
                                 coords=tuple(float(i) for i in line[1:4])))
        return coords

    def write_xyz(self, atoms, filename=None):
        """
        Writes an xyz file using a list of |Atom| instances
        """
        if filename is None:
            raise ValueError('write_xyz: Must give a path to the output file')
        else:
            with open(filename, "w") as f:
                f.write(str(len(atoms)) + '\n\n')
                for atom in atoms:
                    f.write(
                        f"{atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f} \n"
                    )

    def check_db(self):
        """
        Checks fragments for a match in the database
        """

        symbols = {}
        for frag, atoms in self.mol_dict.items():
            symbols[frag] = [atom.symbol for atom in atoms]

        def check_dict(molecules_dict, symbols_dict, db):
            """
            Checking molecule database for a match, and returning the required attributes- 
            name, atoms in molecule, type of molecule, charge, multiplicity, 
            and elements present in the molecule"""
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
            elif db == Molecule.Dications:
                charge = 2
                mult = 1
                mol_type = 'dication'
            elif db == Molecule.Anion_radicals:
                charge = -1
                mult = 2
                mol_type = 'anion-radical'
            elif db == Molecule.Cation_radicals:
                charge = 1
                mult = 2
                mol_type = 'cation-radical'
            elif db == Molecule.Dication_radicals:
                charge = 2
                mult = 2
                mol_type = 'dication-radical'

            for name, atom_list in db.items():
                for sym, molecule in symbols_dict.items():
                    if sorted(molecule) == sorted(atom_list):
                        data = {
                            "type": mol_type,
                            "name": name,
                            "atoms": self.mol_dict[sym],
                            "charge": charge,
                            "multiplicity": mult,
                            "elements": sort_elements(self.mol_dict[sym]),
                            "frag_type": "frag"
                        }
                        molecules_dict[sym] = data
            return molecules_dict

        self.fragments = {}
        for db in (Molecule.Cations, Molecule.Anions, Molecule.Neutrals,
                   Molecule.Radicals, Molecule.Anion_radicals,
                   Molecule.Cation_radicals, Molecule.Dications,
                   Molecule.Dication_radicals):
            self.fragments = check_dict(self.fragments, symbols, db)

        #sort order of atoms
        for data in self.fragments.values():
            data['atoms'] = sorted(data['atoms'], key=lambda atom: atom.index)
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
        convert_keys = {k: v
                        for k, v in enumerate(self.fragments.keys(), 1)
                        }  # old: new
        frags = list(self.fragments.items())
        self.fragments.clear()
        for k, v in frags:
            for key, val in convert_keys.items():
                if k == val:
                    self.fragments[key] = v
                    for atom in self.fragments[key]['atoms']:
                        atom.mol = key

    def give_atoms_a_fragment_name(self):
        for num, frag in self.fragments.items():
            for atom in frag['atoms']:
                atom.fragment = f"{frag['name']}_{num}"

    def sort_fragments_by_index(self):
        """
        Sorts self.fragments according to the index of the atoms in each fragment.
        Sorts by the index of the first atom in each fragment
        """
        frags = [(k, v) for k, v in self.fragments.items()]
        frags = sorted(frags, key=lambda kv: kv[1]['atoms'][0].index)
        self.fragments = {k: v for k, v in frags}

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
                    start, end = 0, 0  # start looking again
                elif start != 0 and end == 0:
                    # found first atom
                    pass
                else:
                    frag_indices.append(int(frag))  # single number

            return frag_indices
            # [(start of frag, end of frag), single atom, (...), (...)]

        def update_mol_dictionary(frag_indices):
            for molecule, item in enumerate(frag_indices):
                num = molecule + 1
                if isinstance(item, tuple):  # if molecule (1, 21) --> choline
                    start, end = item
                    for ind in range(start, end + 1):
                        self.coords[ind - 1].index = ind
                        self.coords[ind - 1].mol = num
                else:  # single-atom 'molecule' i.e. Li, Na, K, Cl, Br
                    self.coords[item - 1].index = item
                    self.coords[item - 1].mol = num

            # reassign mol dict, check db
            self.mol_dict.clear()
            mols = set([atom.mol for atom in self.coords])
            for mol in mols:
                self.mol_dict[mol] = Molecule(
                    atoms=[atom for atom in self.coords if atom.mol == mol])
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
            ionic_mol = Molecule(atoms=coords)
            # causing errors here I think, when separating and using for input
            ionic_mol.separate()
            # instead, just use the mol assignments from the original coords list
            ionic_frags = ionic_mol.fragments
            charge = sum(frag['charge'] for frag in ionic_frags.values())
            multiplicity = Molecule.get_multiplicity(ionic_frags)

            return charge, multiplicity

        coord_list = [coord for coord in self.coords]
        for k, frag in self.fragments.items():
            # remove neutrals
            if frag['charge'] == 0:
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
            # charge, multiplicity = ionic_mol_properties(coord_list)
            self.ionic = {
                "type": 'ionic',
                "name": 'ionic',
                "atoms": coord_list,
                "charge": 0,
                "multiplicity": 1,
                "elements": sort_elements(coord_list),
                "frag_type": "ionic"
            }

    @property
    def all_atoms_assigned(self):
        """
        Checks self.coords to see if all atoms have a molecule assigned
        """
        for atom in self.coords:
            if atom.mol is None:
                return False
        return True

    @property
    def all_fragments_known(self):
        """
        Check that all fragments are in the database
        """

        if not hasattr(self, 'fragments'):
            self.separate()

        for atom in self.coords:
            try:
                name_given = self.fragments[atom.mol]['name']
                if name_given not in Molecules.molecules:
                    return False
            except KeyError:
                # means that fragments haven't been renumbered, so
                # must be something not right
                return False
        return True

        # return all( in molecules for frag in self.fragments.values())

    def group_frags_together(self):
        """
        If the |Molecule| instance has an attribute `group_together`, i.e.
        self.group_together = 'lithium-saccharinate',
        then fragments of lithium and saccharinate are combined into one.
        If more than one of each molecule exists, only one is combined into
        a new fragment at a time- the idea is to combine fragments to run
        jobs on multiple nodes, so that a node doesn't have just one atom
        assigned to it.
        """
        def merged_type(charge):
            if charge == 1:
                return 'cation'
            elif charge == 0:
                return 'neutral'
            elif charge == -1:
                return 'anion'
            else:
                return f'charge: {charge}'

        groups = self.group_together.split('-')
        frags = [frag['name'] for frag in self.fragments.values()]
        merge_keys = []
        for f in frags:
            for k, v in self.fragments.items():
                if v['name'] == f and f in groups:
                    merge_keys.append(k)

        new_key = max(self.fragments.keys()) + 1
        merged_atoms = []
        merged_charge = 0
        merged_elements = []
        for key in merge_keys:
            merged_atoms += self.fragments[key]['atoms']
            merged_charge += self.fragments[key]['charge']
            merged_elements += self.fragments[key]['elements']
        merged_elements = set(merged_elements)
        merged_mult = 2 if any(frag['multiplicity'] == 2
                               for frag in self.fragments.values()) else 1
        merged_mol_type = merged_type(merged_charge)
        self.fragments[new_key] = {
            'type': 'merged',  # was merged_mol_type
            'name': self.group_together,
            'atoms': merged_atoms,
            'charge': merged_charge,
            'multiplicity': merged_mult,
            'elements': merged_elements,
            'frag_type': 'frag'
        }
        for key in merge_keys:
            del self.fragments[key]
        self.fragments_after_merge = self.fragments
        self.frags_grouped_if_desired = True

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
        if not self.all_atoms_assigned:  #or not self.all_fragments_known: # fix for stampede check_hf_v_mp2 geodesics
            print(
                f"{len(self.fragments)} fragments found, some atoms unaccounted for..."
            )
            all_assigned = False
            while not all_assigned:
                self.reassign_frags_manually()
                if self.all_atoms_assigned:
                    all_assigned = True
        if hasattr(self,
                   'group_together') and not self.frags_grouped_if_desired:
            self.group_frags_together()
        self.give_atoms_a_fragment_name()
        self.add_ionic_network()
        if hasattr(self, 'fragments_after_merge'):
            self.fragments = self.fragments_after_merge
        if self.split_on_bonds:
            self.fragment_on_bonds()

    def fragment_on_bonds(self):
        """
        Takes a system that has already been fragmented according to intermolecular
        distance, and then fragments again according to the bonds passed in by the 
        `bonds_to_split` parameter. This should be a nested list of atom indices,
        indicating which bond to break. For example, [(4,9)] indicates a bond between
        atoms 4 and 9 of the original xyz file that should be broken. 
        """
        def remove_connection(connections, atom1, atom2):
            if atom2 in connections[atom1]:
                connections[atom1].remove(atom2)
            if atom1 in connections[atom2]:
                connections[atom2].remove(atom1)
            return connections

        connections = {
            atom.index: [con.index for con in atom.connected_atoms]
            for atom in self
        }
        # apply split
        for bond in self.bonds_to_split:
            a1, a2 = bond
            connections = remove_connection(connections, a1, a2)

        # convert indices back to atom objects
        connections = {
            k: [self.coords[i - 1] for i in v]
            for k, v in connections.items()
        }

        # now have {original_atom: [connections_to_original_atom]}
        # but need include original_atom in that dict

        for k, v in connections.items():
            # include as first, why not
            v.insert(0, self.coords[k - 1])

        # collect up frags after applying split
        for index, connected_atoms in connections.items():
            # now check rest of indices of connections
            # for connections to any of current atoms in fragment
            for index2, connected2 in connections.items():
                if index != index2:
                    if any(atom in connected2
                           for atom in connected_atoms):  # frags are connected
                        # add the rest of connected2 into connected_atoms and empty 'old' frag
                        for a2 in connected2:
                            if a2 not in connected_atoms:
                                connected_atoms.append(a2)
                        connections[index2] = []

        connections = {
            k: sorted(v, key=lambda atom: atom.index)
            for k, v in connections.items() if len(v) != 0
        }

        # redefine molecule number for each atom, starting from 1
        redefined = {}

        for new_index, kv in enumerate(connections.items(), 1):
            key, value = kv
            redefined[new_index] = value

        for k, v in redefined.items():
            for atom in v:
                atom.mol = k

        self.split_fragments = redefined

        # how do we assign charges??????
        # need to look back at original fragments...

        # assume that breaking does not change charge at all, and just carry over the original charges
        # - if original molecule not charged, then nothing to do, set charge to 0, mult to 1
        # - if molecule was originally charged, now can do one of two things:
        #   1) Define certain functional groups to have charges
        #       i.e. SO₃⁻ = -1, RNH3⁺ = +1 etc...
        #   2) Count all electrons of fragment and determine electronic charge from there?
        # For ease of use, can define fragments and read in from file

        # At end of this, reassign self.fragments to be self.split_fragments, so that GAMESS
        # script can use this new fragment dict for the FMO calcs.

        # molecules already read from ~/.config/autochem/molecules.txt
        # so can check if atoms are in the Molecules dict

        # still not perfect because users have to input their fragments manually in the db

        def check_frag_in_db(atoms):
            """
            Checking molecule database for a match, and returns charge and mult 
            in that order.
            If not found, returns a neutral species with no unpaired electrons.
            """
            anions = [sorted(v) for v in Molecule.Anions]
            cations = [sorted(v) for v in Molecule.Cations]
            neutrals = [sorted(v) for v in Molecule.Neutrals]
            rads = [sorted(v) for v in Molecule.Radicals]
            dicats = [sorted(v) for v in Molecule.Dications]
            an_rads = [sorted(v) for v in Molecule.Anion_radicals]
            cat_rads = [sorted(v) for v in Molecule.Cation_radicals]
            dicat_rads = [sorted(v) for v in Molecule.Dication_radicals]
            if atoms in anions:
                return -1, 1
            elif atoms in cations:
                return 1, 1
            elif atoms in neutrals:
                return 0, 1
            elif atoms in rads:
                return 0, 2
            elif atoms in dicats:
                return 2, 1
            elif atoms in an_rads:
                return -1, 2
            elif atoms in cat_rads:
                return 1, 2
            elif atoms in dicat_rads:
                return 2, 2
            else:
                return 0, 1

        #sort order of atoms
        for data in self.fragments.values():
            data['atoms'] = sorted(data['atoms'], key=lambda atom: atom.index)
            # add number
            for i, atom in enumerate(data['atoms']):
                atom.number = i + 1
        # reassigning self.fragments
        self.fragments = {}
        for k, v in self.split_fragments.items():
            atoms = sorted(v, key=lambda atom: atom.index)
            charge, mult = check_frag_in_db(atoms)
            self.fragments[k] = {
                'type': 'frag',
                'name': f'fragmented_{k}',
                'atoms': v,
                'charge': charge,
                'multiplicity': mult,
                'elements': sort_elements(v),
                'frag_type': 'fragmented_on_bond'
            }

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
                    if dists[i, j] < vdw_dist:
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
        self.mol_dict = {
            val: [atom for atom in self.coords if atom.mol == val]
            for val in nums
        }
        for mol in self.mol_dict.values():
            mol.sort(key=lambda atom: atom.index)

    def find_h_bonds(self, distance=2.0):
        """
        Gives hydrogen bonding data back to the user. Checks for suitable
        connected atoms and bond lengths of less than 2 Å, and bond angles of
        45° either side of linear. Users can also provide an optional distance.
        """
        def find_bonds(self):
            def is_imid_c2_h(atom):
                """
                Checks if a C2-H proton of imidazolium is found
                """
                connectors = {}
                if atom.symbol == 'H':
                    for alpha in atom.connected_atoms:  # alpha = one atom away, beta = two away
                        for beta in alpha.connected_atoms:
                            if beta.symbol not in connectors:
                                connectors[beta.symbol] = 1
                            else:
                                connectors[beta.symbol] += 1
                            if alpha.symbol == 'C' and 'N' in connectors and connectors[
                                    'N'] == 2:
                                return True
                return False

            def is_alkyl(atom):
                if atom.symbol == 'H':
                    for a in atom.connected_atoms:
                        if a.symbol == 'C':
                            return True
                return False

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

                # can't have two hydrogens for example
                # if all(atom.symbol is 'H' for atom in (atom1, atom2)):
                if atom1.symbol == atom2.symbol:
                    return False
                for atom in (atom1, atom2):
                    if atom.symbol == 'H':
                        if is_imid_c2_h(atom):  # exception
                            return True
                        if is_alkyl(atom):
                            return False
                    if atom.symbol not in h_bonders:
                        return False
                        for a in atom.connected_atoms:
                            if a.symbol not in h_bonders:
                                return False

                return True

            def within_hbond_distance(atom1, atom2, dist):
                """
                Checks that atoms are within hydrogen-bonding distances, set to
                2 Å by default.
                """
                return atom1.distance_to(atom2) < dist

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
                return atom2.angle_between(atom1.coords,
                                           connected_to_atom2.coords)

            def within_angle_tolerance(atom1, atom2):
                """
                Checks that a hydrogen bond is formed at a suitable angle, 45°
                either side of linear.
                """
                return 225 > bond_angle(atom1, atom2) > 145

            def valid_bond(atom1, atom2, dist):
                """
                Checks that two atoms forms a valid hydrogen bond.
                """
                if within_hbond_distance(atom1, atom2, dist) and \
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
                                if valid_bond(atom1, atom2, distance):
                                    pair = sorted([atom1.index, atom2.index])
                                    if pair not in counted:
                                        dist = atom1.distance_to(atom2)
                                        angle = bond_angle(atom1, atom2)
                                        h_bonded.append(
                                            [atom1, atom2, dist, angle])
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
            db['An'] = Molecule.Anions
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
                # dmso_1
                mol_one_name = f"{self.fragments[one.mol]['name']}_{one.mol}"
                mol_two_name = f"{self.fragments[two.mol]['name']}_{two.mol}"

                group1 = find_molecule_type(mol_one_name)
                group2 = find_molecule_type(mol_two_name)

                hbond_data.append([
                    mol_one_name, one.symbol, mol_two_name, two.symbol, dist,
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
        return 2 if any('radical' in frag['type']
                        for frag in fragment_dict.values()) else 1
        # extend multiplicity for biradicals etc...

    def mol_template(self):
        lines = [
            "# Molecules should be laid out in four lines as follows:\n",
            "# name=<NAME>\n", "# charge=<CHARGE>\n",
            "# multiplicity=<MULTIPLICITY>\n",
            "# atoms=<list of individual atoms in any order>\n",
            "# hashed and blank lines are not read by python,\n",
            '# and names should contain no spaces\n'
            "# below is an example for hydrogen peroxide:\n\n", "name=h2o2\n",
            "charge=0\n", "multiplicity=1\n", "atoms=O,H,H,O\n"
        ]
        return lines

    def check_user_additions(self):
        """
        Reads ~/.config/autochem/molecules.txt for 
        any additional molecules
        """
        confdir = os.path.expanduser('~/.config/autochem/')
        if not os.path.exists(confdir):
            os.makedirs(confdir, exist_ok=True)
        userfile = os.path.join(confdir, 'molecules.txt')
        # CREATE TEMPLATE FILE IF NOT EXISTS
        if not os.path.isfile(userfile):
            open(userfile, 'w+').writelines(self.mol_template())
        # READ USER MOLECULES IF FILE EXISTS
        else:
            name = False
            charge = False
            mult = False
            atoms = False
            with open(userfile, 'r+') as f:
                for line in f:
                    # GET RID OF EXTRA SPACES AND ANYTHING AFTER A HASH
                    line = line.strip()
                    line = line.split('#')[0]
                    # SPLIT INTO DESCRIPTOR AND VALUE
                    line = line.split('=')
                    # FIND IF ONE OF THE DESCRIPTORS AND ASSIGN VALUE
                    if 'name' in line[0]:
                        name = line[1]
                    elif 'charge' in line[0]:
                        charge = int(line[1])
                    elif 'multiplicity' in line[0]:
                        mult = int(line[1])
                    elif 'atoms' in line[0]:
                        atoms = line[1].split(',')
                        for i in range(len(atoms)):
                            atoms[i] = atoms[i].strip()

                    # ONCE ALL DEFINED ADD TO DICTIONARY
                    if not name is False and not charge is False:
                        if not mult is False and not atoms is False:
                            if charge == 0 and mult == 1:
                                Molecule.Neutrals[name] = atoms
                            elif charge == -1 and mult == 1:
                                Molecule.Anions[name] = atoms
                            if charge == 1 and mult == 1:
                                Molecule.Cations[name] = atoms
                            if charge == 0 and mult == 2:
                                Molecule.Radicals[name] = atoms
                            if charge == -1 and mult == 2:
                                Molecule.Anion_radicals[name] = atoms
                            if charge == 1 and mult == 2:
                                Molecule.Cation_radicals[name] = atoms
                            if charge == 2 and mult == 1:
                                Molecule.Dications[name] = atoms
                            if charge == 2 and mult == 2:
                                Molecule.Dication_radicals[name] = atoms
                            # RESET VARS
                            name = False
                            charge = False
                            mult = False
                            atoms = False
