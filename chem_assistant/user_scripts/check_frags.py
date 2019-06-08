__all__ = ['print_frags']

from ..core.molecule import Molecule
import os

def print_frags(directory):
    """
    Prints fragments for each xyz file in the directory passed in
    """
    files = [file for file in os.listdir(directory) if file.endswith('xyz')]

    print('-' * 30)
    for file in files:
        print(file)
        mol = Molecule(using = file)
        mol.separate()
        for frag in mol.fragments.values():
            print(frag['name'], end = '')
            print(f"\t {frag['atoms'][0].index}-{frag['atoms'][-1].index}")
        print('-' * 30)
