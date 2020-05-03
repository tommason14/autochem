from ..core.molecule import Molecule
from ..core.utils import responsive_table
import os

__all__ = ["print_frags"]


def print_frags(directory, verbose=False, grouping=None):
    """
    Prints fragments for each xyz file in the directory passed in.
    If the verbose setting is passed, print out indices of each fragment.
    If grouping is not None, then fragments of those names are grouped 
    together.
    """
    files = [file for file in os.listdir(directory) if file.endswith("xyz")]

    print()
    for file in files:
        mol = Molecule(using=file, group=grouping)
        mol.separate()
        if len(mol.fragments) == 0:
            print(f'{file}: contains no molecules found in the database.\nConsider adding molecules to the ~/.config/autochem/molecules.txt file.\n')
            continue
        else:
            print(f"{file.replace('.xyz', '')}: {len(mol.fragments)} fragments")
            for_printing = {"Fragments": [f["name"] for f in mol.fragments.values()]}
            if verbose:
                for_printing["Atoms"] = [
                    f"{' '.join([str(atom.index) for atom in frag['atoms']])}"
                    for frag in mol.fragments.values()
                ]
            else:
                for_printing["Atoms"] = [
                    f"{frag['atoms'][0].index}-{frag['atoms'][-1].index}"
                    for frag in mol.fragments.values()
                ]

            # for frag in mol.fragments.values():
            #     print(frag["name"], end="\t")
            #     if verbose:
            #         print(f"{' '.join([str(atom.index) for atom in frag['atoms']])}")
            #     else:
            #         print(f"{frag['atoms'][0].index}-{frag['atoms'][-1].index}")
            responsive_table(for_printing, strings=[1, 2])
            print()
