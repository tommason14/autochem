from periodic_table import PeriodicTable as PT
from atom import Atom
from bond import Bond
import numpy as np

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
    Cations['choline'] = ['N', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'O', 'H']

    Neutrals = {"nh3" : ['N', 'H', 'H', 'H']}
    Neutrals['water'] = ['H', 'H', 'O']


    def __init__(self, using):
        self.coords = self.read_xyz(using)
        #list of Atom objects, more useful than list of coordinates


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
        return f"{mass:.2f} g/mol\u207b\u00b9"


    def read_xyz(self, using):
        """Reads coordinates of an xyz file and return a list of |Atom| objects, one for each atom"""
        coords = []
        with open(using, "r") as f:
            for coord in f.readlines()[2:]:
                line = coord.split()
                for val in PT.ptable.values():
                    if line[0] == val[0]:
                        coords.append(Atom(line[0], coords = tuple([float(i) for i in line[1:4]])))
        return coords

    def write_xyz(self, atoms, filename = None):
        """Writes an xyz file using a list of |Atom| instances.
        Params:

        * ``atoms`` -- list of |Atom| objects
        """
        if filename is None:
            raise ValueError('write_xyz: Must give a path to the output file')
        else:
            with open(filename, "w") as f:
                f.write(str(len(atoms)) + '\n\n')
                for atom in atoms:
                    f.write(f"{atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f} \n")



    def separate(self, cutoff = 1.7):
        """Assigns each atom in ``self.coords`` to a different fragment- very useful for FMO calculations"""
        
        def split(self, cutoff):
            # find distance between all atoms, store in a square matrix
            distances = np.zeros((len(self.coords), len(self.coords)))
            group = 0
            for i, atom_i in enumerate(self.coords):
                connected = False
                for j, atom_j in enumerate(self.coords):
                    if atom_i != atom_j:
                        distances[i][j] = atom_i.distance_to(atom_j)
                        # group atoms together based on distance
                        if distances[i, j] < cutoff:
                            connected = True
                            # check if connected using vdw radii
                            if distances[i, j] < PT.get_radius(atom_i.symbol) + PT.get_radius(atom_i.symbol):
                                if atom_j not in atom_i.bonded_atoms:
                                    atom_i.bonded_atoms.append(atom_j)
                                    atom_j.bonded_atoms.append(atom_i)
                                    atom_i.bonds.append(Bond(self, atom_j))
                                    atom_j.bonds.append(Bond(self, atom_i))

                            # IF NEITHER I NOR J PART OF A GROUP
                            # ADD THEM TO A NEW GROUP
                            if atom_i.mol is None and atom_j.mol is None:
                                atom_i.mol, atom_j.mol = group, group
                                group += 1
                            # IF BOTH HAVE BEEN ASSIGNED TO A DIFF GROUP
                            elif not atom_i.mol is None and not atom_j.mol is None:
                                if atom_i.mol != atom_j.mol:
                                    grp_chng = atom_j.mol
                                    for a in self.coords:
                                        if a.mol is grp_chng:
                                            a.mol = atom_i.mol
                            # IF j NOT ASSIGNED
                            elif not atom_i.mol is None and atom_j.mol is None:
                                atom_j.mol = atom_i.mol
                            # IF i NOT ASSIGNED
                            elif not atom_j.mol is False and atom_i.mol is False:
                                atom_i.mol = atom_j.mol
                if not connected:
                    atom_i.mol = group
                    group += 1

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
            """Returns a dictionary of fragments, complete with names, multiplicities and charges, ready for output either to inidividual files, or as input for a GAMESS FMO calculation"""
            molecules = {}

            symbols = {}
            for frag in frags:
                for atom in frags[frag]:
                    if frag not in symbols:
                        symbols[frag] = [atom.symbol]
                    else:
                        symbols[frag].append(atom.symbol)

            def check_dict(molecules, symbols, frags, db, charge, mult, mol_type):
                """Checking molecule database for a match, and returning the required attributes- name, atoms in molecule, type of molecule, charge, multiplicity"""
                for name, atom_list in db.items():
                  for sym, molecule in symbols.items():
                    if sorted(molecule) == sorted(atom_list): #sorts in place, but no overwriting of variable
                        molecules[sym] = {
                            "type": mol_type,
                            "name" : name,
                            "atoms": frags[sym], #not symbols, but the atom objects
                            "charge": charge,
                            "multiplicity": mult
                        }
                return molecules

            check_dict(molecules, symbols, frags, Molecule.Anions, -1, 0, 'anion')
            check_dict(molecules, symbols, frags, Molecule.Cations, 1, 0, 'cation')
            check_dict(molecules, symbols, frags, Molecule.Neutrals, 0, 0, 'neutral')

            return molecules


        def attempt(cutoff = 1.7):
            split(self, cutoff)
            frags = grouping(self)
            molecules = check_db(self, frags)
            return molecules


        def print_fragments(molecules):
            for mol, value in molecules.items():
                print(f"{value['type'].capitalize()} found: {value['name']}")
            correct = input(f"Does your system have {len(molecules)} fragments? [Y/N] ").lower()
            while correct not in ('y', 'n'):
                print("Please type 'Y' or 'N'")
                correct = input(f"Does your system have {len(molecules)} fragments? [Y/N] ").lower()

            return correct

        correct = 'n'
        cutoff = 1.7
        while correct != 'y':
            split(self, cutoff)
            frags = grouping(self)
            molecules = check_db(self, frags)
            correct = print_fragments(molecules)
            if correct != 'y':
                cutoff = float(input('Cutoff distance (Å): '))