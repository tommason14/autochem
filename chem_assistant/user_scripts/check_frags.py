from ..core.molecule import Molecule
import os

__all__ = ['print_frags']

def print_frags(directory, verbose = False, grouping = None):
    """
    Prints fragments for each xyz file in the directory passed in.
    If the verbose setting is passed, print out indices of each fragment.
    If grouping is not None, then fragments of those names are grouped 
    together.
    """
    files = [file for file in os.listdir(directory) if file.endswith('xyz')]

    print('-' * 30)
    for file in files:
        print(file)
        mol = Molecule(using=file, group = grouping)
        mol.separate()
        for frag in mol.fragments.values():
            print(frag['name'], end='')
            if verbose:
                print(f"\n\t {' '.join([str(atom.index) for atom in frag['atoms']])}")
            else:
                print(f"\t {frag['atoms'][0].index}-{frag['atoms'][-1].index}")
        print('-' * 30)
