#!/usr/bin/env python3

"""
File: allocate_vibrations.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Takes the output from gamess_to_molden.py
and attempts to assign each vibration to the atom in question.
This involves splitting up the system into molecules, 
finding atoms that are bonded, and then matching the atoms that
move the most with the bond that connects the moving atom
to any others.

Note: If you are predicting vibrations of a new system, add your molecule into
the dictionaries starting on line 360.

Code that performs any function starts on line 780

TODO: Ranking likelihood of each vibration, with cli -r number
i.e. -r 3 for the top 3 of each vib
"""
import re
import numpy as np
import math
import itertools
import sys
import argparse

parser = argparse.ArgumentParser(description='Make a guess at vibrations present in molden output files')
# required = parser.add_argument_group('required arguments')
parser.add_argument('-f', '--file', help='Molden file to use', action='store', required = True)
parser.add_argument('-i', '--sort-by-intensity', help='Sorts output by intensity of vibrations', action='store_true')
parser.add_argument('-l', '--low', help='Lowest frequency to look for. Pass in an integer', action='store', type = int)
parser.add_argument('-m', '--high', help='Highest frequency to look for (m for max). Pass in an integer', action='store', type = int)
parser.add_argument('-w', '--weak', help='Lowest intensity transition to look for (w for weak). Pass in an integer', action='store', type = int)
parser.add_argument('-s', '--strong', help='Highest intensity transition to look for (s for strong). Pass in an integer', action='store', type = int)
parser.add_argument('-n', '--number', help='Print the first n vibrations', action='store', type = int)
parser.add_argument('-t', '--tail', help='Print the last t vibrations', action='store', type = int)
args = parser.parse_args()

def responsive_table(data, strings, min_width):
    """
    Returns a table that is reponsive in size to every column.
    Requires a dictionary to be passed in, with the keys referring to
    the headers of the table.
    Also pass in the number of each column that should be a string, starting
    from 1.

    Usage:
        >>> d = {'col1': [1,2,3,4],
                 'col2': ['One', 'Two', 'Three', 'Four']}
        >>> responsive_table(d, strings = [2])

    Can also give a minimum width, defaults to 13 spaces
    """
    num_cols = len(data.keys())
    content = zip(*[data[key] for key in data.keys()]) # dict values into list of lists
    # unknown number of arguments
    max_sizes = {}
    for k, v in data.items():
        max_sizes[k] = len(max([str(val) for val in v], key = len))

    # create the thing to pass into .format()- can't have brackets like zip gives
    formatting = []
    index = 0
    all_sizes = []
    if min_width is None:
        min_width = 13
    for val in zip(data.keys(), max_sizes.values()):
        entry, size = val
        if size < min_width or index + 1 not in strings:
            size = min_width
        # also check dict key length
        if len(entry) > size:
            size = len(entry)
        formatting += [entry, size]
        all_sizes.append(size)
        index += 1
    line_length = sum(all_sizes) + num_cols * 3 - 1 # spaces in header
    print('+' + '-' * line_length + '+')
    output_string = '|' + " {:^{}} |" * len(data.keys())
    print(output_string.format(*formatting))
    print('+' + '-' * line_length + '+')
    for line in content:
        formatting = []
        for val in zip(line, all_sizes):
            entry, size = val
            if not isinstance(entry, str):
                size = f'{size}.2f'
            formatting.append(entry)
            formatting.append(size)
        print(output_string.format(*formatting))
    print('+' + '-' * line_length + '+')


def read_file(file):
    with open(file, "r") as f:
        try:
            for line in f:
                yield line
        except UnicodeDecodeError:
            pass

def sort_elements(lst):
    """
    Sort a list of |Atom| objects by atomic number.
    Returns a list of tuples- [(symbol, atomic number), (symbol, atomic number), ...]

    TODO: Extend to giving back the objects- more useful than just for formatting of symbols
    """
    els = []
    elements = set([atom.symbol for atom in lst])
    for i in elements:
        atom = Atom(i)
        els.append((i, float(PT.get_atnum(atom))))
    sorted_els = sorted(els, key = lambda val: val[1])
    return sorted_els


class Atom:
    """A class representing a atom in 3 dimensional euclidean space.
    An instance has the following attributes:

    * ``atnum`` -- atomic symbol, equal to zero for a dummy atom
    * ``coords`` -- tuple of x,y,z coordinates
    * ``bonds`` -- list of bonds that this atom is a part of
    * ``mol`` -- molecule this atom is a part of. Assigned programmatically when a molecule is separated using the *mol.separate* method, or can be assigned manually if building up a molecule from scratch

    Access these properties directly:
    * ``x``, ``y``, ``z`` -- for atom coordinates
    * ``symbol`` -- read or write the atom symbol directly

    >>> a = Atom('H', coords = (1,2,3))

    """
    def __init__(self, symbol = None, atnum = 0, coords = None, mol = None, bonds = None):
        if symbol is not None:
            self.symbol = symbol
            self.atnum = PT.get_atnum(self)
        else:
            self.atnum = atnum
            self.symbol = PT.get_symbol(self)

        self.mol = mol
        self.bonds = bonds or []
        self.connected_atoms = []
        self.h_bonded_to = []

        if coords is None:
            self.coords = (0, 0, 0)
        elif len(coords) == 3:
            self.coords = tuple(float(i) for i in coords)
        else:
            raise TypeError('Atom: Invalid coordinates given')
        self.x, self.y, self.z = self.coords

    def __repr__(self):
        """Unambiguous representation of an |Atom| instance"""
        if hasattr(self, 'index') and not hasattr(self, 'mol'):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Index: {self.index} "
        if hasattr(self, 'index') and hasattr(self, 'mol'):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Index: {self.index} Mol: {self.mol}"
        elif hasattr(self, 'number'):
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Mol: {self.mol} Atom: {self.number}"
        elif hasattr(self, 'index') and len(self.h_bonded_to) > 0:
            h_bonded = [(atom.symbol, {'mol': atom.mol, 'atom': atom.index}) for atom in self.h_bonded_to]
            return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}\
 Mol: {self.mol} Index: {self.index} Number: {self.number} H-Bonds: {h_bonded}"
        return f"Atom: {self.symbol:3s} {self.x:>10.5f} {self.y:>10.5f} {self.z:>10.5f}"

    def __iter__(self):
        """Iterates through coordinates when called"""
        return iter(self.coords)

    def translate(self, vector):
        """Move atom in space by passing a vector in angstroms"""
        self.coords = tuple(i + j for i,j in zip(self, vector))

    def move_to(self, vector):
        """Move atom in space to the values, in angstroms, given in this vector. The vector passed represents a point in euclidean space"""
        self.coords = tuple(i for i in vector)

    def distance_to(self, vector):
        """Measure the distance between the atom and a point in space, given as a vector in angstroms"""
        # pythagoras in 3D
        dist = 0.0
        for i,j in zip(self, vector):
            dist += (i - j)**2
        return dist ** 0.5

    def vector_to(self, point):
        """Returns a vector from the atom to a given point, in angstroms"""
        return tuple((i - j) for i,j in zip(point, self))

    def angle_between(self, pos1, pos2):
        """Returns an angle between positions 1 and 2 in degrees, with this atom lying at the centre"""
        # dot product, angle = cos^-1([vec(a).vec(b)] / [dist(a) * dist(b)])
        num = np.dot(self.vector_to(pos1), self.vector_to(pos2))
        denom = self.distance_to(pos1) * self.distance_to(pos2)
        return math.acos(num/denom) * (180 / math.pi)

    __str__ = __repr__


class PT:
    """Helper class to allow for lookup of atomic properties. Can convert between symbol and atomic number"""
    ptable = {}
    #[symbol, mass, radius, connectors, vdw radii]
        #atomic weights from: http://www.ciaaw.org/atomic-weights.htm
    ptable[  0] = ['Xx',   0.00000, 0.00 ,  0, 0.000]
    ptable[  1] = [ 'H',   1.00798, 0.30 ,  1, 0.430]
    ptable[  2] = ['He',   4.00260, 0.99 ,  0, 0.741]
    ptable[  3] = ['Li',   6.96750, 1.52 ,  8, 0.880]
    ptable[  4] = ['Be',   9.01218, 1.12 ,  8, 0.550]
    ptable[  5] = [ 'B',  10.81350, 0.88 ,  6, 1.030]
    ptable[  6] = [ 'C',  12.01060, 0.77 ,  4, 0.900]
    ptable[  7] = [ 'N',  14.00685, 0.70 ,  3, 0.880]
    ptable[  8] = [ 'O',  15.99940, 0.66 ,  2, 0.880]
    ptable[  9] = [ 'F',  18.99840, 0.64 ,  1, 0.840]
    ptable[ 10] = ['Ne',  20.17970, 1.60 ,  0, 0.815]
    ptable[ 11] = ['Na',  22.98977, 1.86 ,  8, 1.170]
    ptable[ 12] = ['Mg',  24.30550, 1.60 ,  8, 1.300]
    ptable[ 13] = ['Al',  26.98154, 1.43 ,  8, 1.550]
    ptable[ 14] = ['Si',  28.08500, 1.17 ,  8, 1.400]
    ptable[ 15] = [ 'P',  30.97376, 1.10 ,  8, 1.250]
    ptable[ 16] = [ 'S',  32.06750, 1.04 ,  2, 1.220]
    ptable[ 17] = ['Cl',  35.45150, 0.99 ,  1, 1.190]
    ptable[ 18] = ['Ar',  39.94800, 1.92 ,  0, 0.995]
    ptable[ 19] = [ 'K',  39.09830, 2.31 ,  8, 1.530]
    ptable[ 20] = ['Ca',  40.07800, 1.97 ,  8, 1.190]
    ptable[ 21] = ['Sc',  44.95591, 1.60 ,  8, 1.640]
    ptable[ 22] = ['Ti',  47.86700, 1.46 ,  8, 1.670]
    ptable[ 23] = [ 'V',  50.94150, 1.31 ,  8, 1.530]
    ptable[ 24] = ['Cr',  51.99610, 1.25 ,  8, 1.550]
    ptable[ 25] = ['Mn',  54.93804, 1.29 ,  8, 1.555]
    ptable[ 26] = ['Fe',  55.84500, 1.26 ,  8, 1.540]
    ptable[ 27] = ['Co',  58.93319, 1.25 ,  8, 1.530]
    ptable[ 28] = ['Ni',  58.69340, 1.24 ,  8, 1.700]
    ptable[ 29] = ['Cu',  63.54600, 1.28 ,  8, 1.720]
    ptable[ 30] = ['Zn',  65.38000, 1.33 ,  8, 1.650]
    ptable[ 31] = ['Ga',  69.72300, 1.41 ,  8, 1.420]
    ptable[ 32] = ['Ge',  72.63000, 1.22 ,  8, 1.370]
    ptable[ 33] = ['As',  74.92159, 1.21 ,  8, 1.410]
    ptable[ 34] = ['Se',  78.97100, 1.17 ,  8, 1.420]
    ptable[ 35] = ['Br',  79.90400, 1.14 ,  1, 1.410]
    ptable[ 36] = ['Kr',  83.79800, 1.97 ,  0, 1.069]
    ptable[ 37] = ['Rb',  85.46780, 2.44 ,  8, 1.670]
    ptable[ 38] = ['Sr',  87.62000, 2.15 ,  8, 1.320]
    ptable[ 39] = [ 'Y',  88.90584, 1.80 ,  8, 1.980]
    ptable[ 40] = ['Zr',  91.22400, 1.57 ,  8, 1.760]
    ptable[ 41] = ['Nb',  92.90637, 1.41 ,  8, 1.680]
    ptable[ 42] = ['Mo',  95.95000, 1.36 ,  8, 1.670]
    ptable[ 43] = ['Tc',  98.00000, 1.35 ,  8, 1.550]
    ptable[ 44] = ['Ru', 101.07000, 1.33 ,  8, 1.600]
    ptable[ 45] = ['Rh', 102.90550, 1.34 ,  8, 1.650]
    ptable[ 46] = ['Pd', 106.42000, 1.38 ,  8, 1.700]
    ptable[ 47] = ['Ag', 107.86820, 1.44 ,  8, 1.790]
    ptable[ 48] = ['Cd', 112.41400, 1.49 ,  8, 1.890]
    ptable[ 49] = ['In', 114.81800, 1.66 ,  8, 1.830]
    ptable[ 50] = ['Sn', 118.71000, 1.62 ,  8, 1.660]
    ptable[ 51] = ['Sb', 121.76000, 1.41 ,  8, 0.000]
    ptable[ 52] = ['Te', 127.60000, 1.37 ,  8, 0.000]
    ptable[ 53] = [ 'I', 126.90447, 1.33 ,  1, 0.000]
    ptable[ 54] = ['Xe', 131.29300, 2.17 ,  0, 0.000]
    ptable[ 55] = ['Cs', 132.90545, 2.62 ,  8, 0.000]
    ptable[ 56] = ['Ba', 137.32700, 2.17 ,  8, 0.000]
    ptable[ 57] = ['La', 138.90547, 1.88 ,  8, 0.000]
    ptable[ 58] = ['Ce', 140.11600, 1.818,  8, 0.000]
    ptable[ 59] = ['Pr', 140.90766, 1.824,  8, 0.000]
    ptable[ 60] = ['Nd', 144.24200, 1.814,  8, 0.000]
    ptable[ 61] = ['Pm', 145.00000, 1.834,  8, 0.000]
    ptable[ 62] = ['Sm', 150.36000, 1.804,  8, 0.000]
    ptable[ 63] = ['Eu', 151.96400, 2.084,  8, 0.000]
    ptable[ 64] = ['Gd', 157.25000, 1.804,  8, 0.000]
    ptable[ 65] = ['Tb', 158.92535, 1.773,  8, 0.000]
    ptable[ 66] = ['Dy', 162.50000, 1.781,  8, 0.000]
    ptable[ 67] = ['Ho', 164.93033, 1.762,  8, 0.000]
    ptable[ 68] = ['Er', 167.25900, 1.761,  8, 0.000]
    ptable[ 69] = ['Tm', 168.93422, 1.759,  8, 0.000]
    ptable[ 70] = ['Yb', 173.04500, 1.922,  8, 0.000]
    ptable[ 71] = ['Lu', 174.96680, 1.738,  8, 0.000]
    ptable[ 72] = ['Hf', 178.49000, 1.57 ,  8, 0.000]
    ptable[ 73] = ['Ta', 180.94788, 1.43 ,  8, 0.000]
    ptable[ 74] = [ 'W', 183.84000, 1.37 ,  8, 0.000]
    ptable[ 75] = ['Re', 186.20700, 1.37 ,  8, 0.000]
    ptable[ 76] = ['Os', 190.23000, 1.34 ,  8, 0.000]
    ptable[ 77] = ['Ir', 192.21700, 1.35 ,  8, 0.000]
    ptable[ 78] = ['Pt', 195.08400, 1.38 ,  8, 0.000]
    ptable[ 79] = ['Au', 196.96657, 1.44 ,  8, 0.000]
    ptable[ 80] = ['Hg', 200.59200, 1.52 ,  8, 0.000]
    ptable[ 81] = ['Tl', 204.38350, 1.71 ,  8, 0.000]
    ptable[ 82] = ['Pb', 207.20000, 1.75 ,  8, 0.000]
    ptable[ 83] = ['Bi', 208.98040, 1.70 ,  8, 0.000]
    ptable[ 84] = ['Po', 209.00000, 1.40 ,  8, 0.000]
    ptable[ 85] = ['At', 210.00000, 1.40 ,  1, 0.000]
    ptable[ 86] = ['Rn', 222.00000, 2.40 ,  0, 0.000]
    ptable[ 87] = ['Fr', 223.00000, 2.70 ,  8, 0.000]
    ptable[ 88] = ['Ra', 226.00000, 2.20 ,  8, 0.000]
    ptable[ 89] = ['Ac', 227.00000, 2.00 ,  8, 0.000]
    ptable[ 90] = ['Th', 232.03770, 1.79 ,  8, 0.000]
    ptable[ 91] = ['Pa', 231.03588, 1.63 ,  8, 0.000]
    ptable[ 92] = [ 'U', 238.02891, 1.56 ,  8, 0.000]
    ptable[ 93] = ['Np', 237.00000, 1.55 ,  8, 0.000]
    ptable[ 94] = ['Pu', 244.00000, 1.59 ,  8, 0.000]
    ptable[ 95] = ['Am', 243.00000, 1.73 ,  8, 0.000]
    ptable[ 96] = ['Cm', 247.00000, 1.74 ,  8, 0.000]
    ptable[ 97] = ['Bk', 247.00000, 1.70 ,  8, 0.000]
    ptable[ 98] = ['Cf', 251.00000, 1.86 ,  8, 0.000]
    ptable[ 99] = ['Es', 252.00000, 1.86 ,  8, 0.000]
    ptable[100] = ['Fm', 257.00000, 2.00 ,  8, 0.000]
    ptable[101] = ['Md', 258.00000, 2.00 ,  8, 0.000]
    ptable[102] = ['No', 259.00000, 2.00 ,  8, 0.000]
    ptable[103] = ['Lr', 266.00000, 2.00 ,  8, 0.000]
    ptable[104] = ['Rf', 267.00000, 2.00 ,  8, 0.000]
    ptable[105] = ['Db', 268.00000, 2.00 ,  8, 0.000]
    ptable[106] = ['Sg', 269.00000, 2.00 ,  8, 0.000]
    ptable[107] = ['Bh', 270.00000, 2.00 ,  8, 0.000]
    ptable[108] = ['Hs', 277.00000, 2.00 ,  8, 0.000]
    ptable[109] = ['Mt', 278.00000, 2.00 ,  8, 0.000]
    ptable[110] = ['Ds', 281.00000, 2.00 ,  8, 0.000]
    ptable[111] = ['Rg', 282.00000, 2.00 ,  8, 0.000]
    ptable[112] = ['Cn', 285.00000, 2.00 ,  8, 0.000]
    ptable[113] = ['Nh', 286.00000, 2.00 ,  8, 0.000]
    ptable[114] = ['Fl', 289.00000, 2.00 ,  8, 0.000]
    ptable[115] = ['Mc', 290.00000, 2.00 ,  8, 0.000]
    ptable[116] = ['Lv', 293.00000, 2.00 ,  8, 0.000]
    ptable[117] = ['Ts', 294.00000, 2.00 ,  8, 0.000]
    ptable[118] = ['Og', 294.00000, 2.00 ,  8, 0.000]
                                               
    def __init__(self): 
        raise AttributeError('The PeriodicTable class cannot be instantiated.')

    @classmethod
    def get_atnum(cls, atom):
        """Converts symbol to atomic number"""
        for key in cls.ptable.keys():
            if cls.ptable[key][0] == atom.symbol.capitalize():
                return key

    @classmethod
    def get_symbol(cls, atom):
        """Converts atomic number to symbol"""
        return cls.ptable[atom.atnum][0]
    
    @classmethod
    def get_radius(cls, atom):
        """Returns atomic radius for a given element"""
        atnum = cls.get_atnum(atom)
        return float(cls.ptable[atom.atnum][2])
    
    @classmethod
    def get_mass(cls, atom):
        """Returns atomic mass for a given element"""
        atnum = cls.get_atnum(atom.symbol)
        return float(cls.ptable[atom.atnum][1])

    @classmethod
    def get_connectors(cls, atom):
        """Returns number of possible attachments to a given element"""
        atnum = cls.get_atnum(atom)
        return int(cls.ptable[atom.atnum][-2])
    
    @classmethod
    def get_vdw(cls, atom):
        """Returns van der waals radius of a given element"""
        atnum = cls.get_atnum(atom)
        return float(cls.ptable[atom.atnum][-1])
        
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

            # list of atoms might not be atom objects
            if not isinstance(atoms[0], Atom):
                new = []
                for atom in atoms:
                        sym, *coords = atom
                        a = Atom(symbol = sym, coords = coords)
                        new.append(a)
                self.coords = new
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
                        if atom_j not in atom_i.connected_atoms:
                            atom_i.connected_atoms.append(atom_j)
                        if atom_i not in atom_j.connected_atoms:
                            atom_j.connected_atoms.append(atom_i)
                    
    def separate(self):
        """
        Separates coordinates into specific fragments using the intermolecular 
        distances along with van der waals radii. Note this function only works 
        with intermolecular fragments and cannot split molecules on bonds.
        """
        self.split_vdw()
        self.check_db()
        self.sort_fragments_by_index()
        self.renumber_molecules()


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

    def frag_name(self, atom):
        """
        Accepts an atom from the molecule, returning the name of the fragment 
        containing that atom.
        """

        if not hasattr(self, 'fragments'):
            self.separate()

        for frag in self.fragments.values():
            for a in frag['atoms']:
                if a.index == atom.index:
                    return frag['name']

#### CODE STARTS HERE

def check_inputs(args):
    """
    Checks if command line arguments are valid
    """
    if len(sys.argv) < 2:
        sys.exit(parser.print_help())
    
    if args.number == args.tail != None: # if defined
        sys.exit("Error: Number of lines from top (n) and number of lines from bottom (t) can't be the same! Change values")

def get_geometry(log):
    """
    Extract the initial geometry from the molden file.

    i.e.
        [Atoms]
            ...
            ... <- take these
            ...
        [FREQ]
    """
    atoms = []
    found = False
    for line in read_file(log):
        if '[Atoms]' in line:
            found = True
        elif '[FREQ]' in line:
            break
        if found:
            line = line.split()
            if len(line) > 1:
                sym = line[0]
                x, y, z = line[-3:]
                atoms.append([sym, x, y, z]) 
    return atoms

def molecule(atoms):
    """
    Creates a molecule object and returns fragments, using
    van der waals radii to split the system
    """
    mol = Molecule(atoms = atoms)
    mol.separate()
    mol.assign_neighbours()
    return mol

def read_freqs(log):
    """
    Obtains frequencies from the molden output file.
    Assumes that we are simulating non-linear polyatomics, returning 
    3N-6 vibrations for N atoms (ignoring symmetry).
    """
    freqs = []
    found = False
    for line in read_file(log):
        if '[FREQ]' in line:
            found = True
        if '[INT]' in line:
            break
        if found:
            try:
                val = float(line.strip()) 
                freqs.append(val)
            except ValueError:
                continue # first line is [FREQ], don't want that
    return freqs[6:]

def read_intensities(log):
    """
    Obtains intensities from the molden output file.
    Assumes that we are simulating non-linear polyatomics, returning 
    3N-6 intensites for N atoms (ignoring symmetry).
    """
    ints = []
    found = False
    for line in read_file(log):
        if '[INT]' in line:
            found = True
        if '[FR-COORD]' in line:
            break
        if found:
            try:
                val = float(line.strip()) 
                ints.append(val)
            except ValueError:
                continue # first line is [INT], don't want that
    return ints[6:]

def read_vibrations(log):
    """
    Obtains the vibrations from the molden output file, in the form of
    displacements from the equilibrium geometry, given as the initial geometry.
    Assumes a non-linear polyatomic molecule, skipping the first 6 vibrations
    and returning 3N-6 in total. 

    Each vibration is returned as a list of x,y,z displacements for each atom,
    with each vibration returned individually. To use, iterate over the
    function as with any iterable:
        >>> for vib in read_vibrations(log):
        >>>     do_something()

    Vibrations are not returned as one list as the number of vibrations could
    become large as the system size increases to many atoms. 
    """
    found = False
    vib_count = 0
    vibs = []
    for line in read_file(log):
        if '[FR-NORM-COORD]' in line:
            found = True
        if found:
            if 'vibration' in line:
                vib_count += 1
            if vib_count > 5:
                if 'vibration' not in line:
                    xyz_list = line.strip().split() 
                    xyz_list = [float(d) for d in xyz_list]
                    vibs.append(xyz_list)
                if 'vibration' in line and len(vibs) > 0: 
                    yield vibs
                    vibs = []

def average_displacement(xyz):
    """
    Accepts a list of x, y, z displacements and returns the average, computing
    the pythagorean distance of the atom from its equilibrium position.
    """
    return sum(disp ** 2 for disp in xyz) ** 0.5

def print_vib(atom_of_max_disp, connected, molecule, sure = True):
    """
    Conventionally, chemists print heteroatoms first when describing bonds
    involving hydrogen. Also takes in a boolean keyword argument of `sure`, for
    if the user is sure that a particular bond is involved, defaulting to True. 
    """

    fragment = molecule.frag_name(atom_of_max_disp)
    bond = ''
    if atom_of_max_disp.symbol == 'H':
        # can make neater, place in dict
        if sure:
            bond = f'{connected.symbol}-{atom_of_max_disp.symbol}'
        else:
            bond = f'Maybe {connected.symbol}-{atom_of_max_disp.symbol}'
    else:
        if sure:
            bond = f'{atom_of_max_disp.symbol}-{connected.symbol}'
        else:
            bond = f'Maybe {atom_of_max_disp.symbol}-{connected.symbol}'
    
    return bond, fragment

def all_same(lst):
    """ 
    Returns true if all elements in the list are the same 
    """
    return len(set(lst)) == 1

def one_outcome(lst):
    """
    Returns true if list only has one item, or if all items are the same
    """
    return len(lst) == 1 or all_same(lst)

def print_vib_info(atom_of_max_disp, molecule):
    """
    Prints the atom(s) that may be involved in the vibration in question, taking
    in the atom of maximum displacement.
    """
    bonds = []
    frags = []

    one_type_of_bond =\
    one_outcome([atom.symbol for atom in atom_of_max_disp.connected_atoms])    

    if one_type_of_bond:
        connected_atom = atom_of_max_disp.connected_atoms[0]
        bond, fragment = print_vib(atom_of_max_disp, connected_atom, molecule)
        bonds.append(bond)
        frags.append(fragment)

    else:
        for connected_atom in atom_of_max_disp.connected_atoms:
            bond, fragment = print_vib(atom_of_max_disp, connected_atom, molecule, sure = False)
            bonds.append(bond)
            frags.append(fragment)
        # now account for displacement direction
        # to find correct bond

    return bonds, frags

def dict_to_nested_list(data):
    """
    Returns a list of lists for a dictionary one level deep. Note dict values
    have to be the same length for this to work properly
    i.e. 
    For a dictionary d,

    d = {'one': [1,2,3,4],
         'two': [5,6,7,8]}

    this function returns
    
    [[1,5], [2,6], [3,7], [4,8]]
    """
    numrows = max(len(v) for v in data.values())
    rows = []
    for num in range(numrows):
        row = []
        for k, v in data.items():
            row.append(v[num])
        rows.append(row)
    return rows

def nested_list_to_dict(lst, keys):
    """
    Takes a nested list and creates a dictionary using the keys attribute, which
    is a list of the desired keys of the dictionary
    """
    data = {}
    for i, value in enumerate(keys):
        data[value] = [row[i] for row in lst]
    return data

def sort_data_by_key(data, keys = None, sort_by=None):
    """
    Sorts a dictionary by second key of `keys`. Sort by the column stord in the
    `sort_by` variable, with the first key being 0. i.e sort_by = 1 sorts on the
    second column
    """
    if sort_by is None:
        raise ValueError('sort_data_by_key: must assign a value to the "sort_by" argument')
    if keys is None:
        keys = data.keys()
    rows = dict_to_nested_list(data)
    rows = sorted(rows, key = lambda row: row[sort_by], reverse = True)
    data = nested_list_to_dict(rows, keys)
    
    return data

def limit_by_key(data, low, high, keys = None, column = None):
    """
    Limit the data stored as values of the dictionary `data`.
    Limit by values `low` and `high`, in the column `column`, where column is
    equal to the position in the key of the dictionary, starting at 0. 
    """
    if keys is None:
        keys = data.keys()
    rows = dict_to_nested_list(data)
    rows = [row for row in rows if low <= row[column] <= high]
    data = nested_list_to_dict(rows, keys)
    return data

def head(data, n):
    """
    Return the first n elements of lists stored as a dictionary value
    """
    keys = data.keys()
    rows = dict_to_nested_list(data)
    rows = rows[:n]
    data = nested_list_to_dict(rows, keys)
    return data

def tail(data, n):
    """
    Return the last n elements of lists stored as a dictionary value
    """
    keys = data.keys()
    rows = dict_to_nested_list(data)
    rows = rows[n:]
    data = nested_list_to_dict(rows, keys)
    return data

def apply_limitations(data, args):

    if args.sort_by_intensity:
        data = sort_data_by_key(data, sort_by=1)

    ### LIMIT FREQUENCIES

    if args.low and not args.high:
        high = max(data['Frequencies (cm⁻¹)'])
        data = limit_by_key(data, args.low, high, column = 0)
    
    if args.high and not args.low:
        low = min(data['Frequencies (cm⁻¹)'])
        data = limit_by_key(data, low, args.high, column = 0)

    if args.low and args.high:
        data = limit_by_key(data, args.low, args.high, column = 0)
    
    ### LIMIT INTENSITIES

    if args.weak and not args.strong:
        strong = max(data['Intensities (au)'])
        data = limit_by_key(data, args.weak, strong, column = 1)

    if args.strong and not args.weak:
        weak = min(data['Intensities (au)'])
        data = limit_by_key(data, weak, args.strong, column = 1)

    if args.strong and args.weak:
        data = limit_by_key(data, args.weak, args.strong, column = 1)

    ### LIMIT ROWS OF OUTPUT
    ### ORDER IMPORTANT IF CUTTING FROM TOP AND BOTTOM

    if args.number and args.tail:
        if args.number > args.tail:
            data = head(data, args.number)
            data = tail(data, args.tail)
        else:
            data = tail(data, args.tail)
            data = head(data, args.number)

    if args.number:
        data = head(data, args.number)

    if args.tail:
        data = tail(data, args.tail)

    return data

def process_vibrations(mol, args):
    """
    Searches through the molden output file for vibrations, then attempts to 
    allocate vibrations to a particular atom by seeing which atom moves the
    most.
    
    Averages displacement across all 3 dimensions to find the atom of largest 
    displacement.
    """

    frequencies = [] 
    intensities = []
    bonds = []
    fragments = []

    freqs = read_freqs(args.file)
    intensity = read_intensities(args.file)

    # more than one vib per frequency, so same frequency appears more than once,
    # need to add same number of intensities
    freq_ints = {}
    for ind, f in enumerate(freqs):
        freq_ints[f] = intensity[ind] 

    for vib, displacements in enumerate(read_vibrations(args.file)):

        freqs_per_vib = []
        bonds_per_vib = []
        frags_per_vib = []
        # print(f'Frequency: {freqs[vib]}')

        all_displacements = {}
        for atom, atom_disp in enumerate(displacements):
            ave_disp = average_displacement(atom_disp)
            all_displacements[atom] = ave_disp

        max_disp = max(all_displacements.values())

        atom_num_of_max_disp = 0
        for k,v in all_displacements.items():
            if v == max_disp:
                atom_num_of_max_disp = k 

        # more than one bond can vibrate per vibration, could just pass in
        # more than one atom of max disp
        # also still need ranking 
        atom_of_max_disp = mol.coords[atom_num_of_max_disp] 
        bonds_max_disp, fragments_max_disp = print_vib_info(atom_of_max_disp, mol)


        bonds_per_vib += bonds_max_disp
        frags_per_vib += fragments_max_disp

        for _ in bonds_per_vib:
            freqs_per_vib.append(freqs[vib])

        frequencies += freqs_per_vib
        bonds += bonds_per_vib
        fragments += frags_per_vib

    # match number of intensities to number of frequencies
    # more frequencies than in log file, if more than one bond is 
    # vibrating
    intensities = []
    for f in frequencies:
        intensities.append(freq_ints[f])
    ### print_table
    
    data = {}
    data['Frequencies (cm⁻¹)'] = frequencies
    data['Intensities (au)'] = intensities
    data['Vibrating bond'] = bonds
    data['Molecule'] = fragments

    data = apply_limitations(data, args)
    responsive_table(data, strings = [2,3], min_width = 10)

def main(args):
    check_inputs(args)
    initial_geom = get_geometry(args.file)
    mol = molecule(initial_geom)
    process_vibrations(mol, args)

if __name__ == '__main__':
    main(args)