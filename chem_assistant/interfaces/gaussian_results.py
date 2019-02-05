from ..core.utils import read_file, write_xyz
from ..core.results import Results

import re
import os
import subprocess

__all__ = ['GaussianResults']

class GaussianResults(Results):
    """Class for obtaining results from Gaussian simulations. This class requires a log file to be read.
    
    MORE TO COME
    """

    def __init__(self, log):
        super().__init__(log)

    ################################
    #                              #
    #      CHECK IF COMPLETED      #
    #                              #
    ################################

    def __repr__(self):
        return self.log # debugging

    __str__ = __repr__

    def get_runtype(self):
        for line in self.read():
            if '#P' in line.upper():
                parts = line.split()
                return parts[2].lower()
                    # if '' in p:
                    #     return p.split('=')[1].lower()

    def completed(self):
        found = False
        for line in self.read():
            if 'Normal termination' in line:
                found = True
        return found

    def get_equil_coords(self, output = None):
        found_equil = False
        found_some = False
        regex = "^(\s*[0-9]*){3}(\s*-?[0-9]{1,3}.[0-9]*){3}$"
        coords = []
        some_coords = []
        par_dir = []
        print(self.log)
        for part in self.path.split('/'):
            if part not in ('opt', 'spec', 'hess'):
                par_dir.append(part)
            else:
                break
        MOLECULE_PARENT_DIR = '/'.join(par_dir)
        for line in self.read():
            if 'Standard orientation' in line:
                found_some = True
            if 'Rotational constants' in line:
                found_some = False 
            if 'Optimization completed' in line:
                found_equil = True
            if 'Distance matrix (angstroms)' in line:
                found_equil = False
            if found_equil:
                if re.search(regex, line):
                    _, atnum, _, x, y, z = line.split()
                    atnum, x, y, z = map(float, (atnum, x, y, z))
                    sym = PT.get_symbol(atnum)
                    coords.append(Atom(sym, coords = (x, y, z)))
            if found_some:
                if len(some_coords) > 0: # from last run, remove those coords
                    some_coords = []
                if re.search(regex, line):
                    _, atnum, _, x, y, z = line.split()
                    atnum, x, y, z = map(float, (atnum, x, y, z))
                    sym = PT.get_symbol(atnum)
                    some_coords.append(Atom(sym, coords = (x, y, z)))
        if len(coords) > 0:
            print('found!')
            write_xyz(coords, os.path.join(MOLECULE_PARENT_DIR, 'equil.xyz'))
        else:
            if len(some_coords) > 0:
                print(f'not found.\nNeeds resubmitting. Coords stored in {self.path}/rerun/rerun.xyz')
                rerun_dir = os.path.join(self.path, 'rerun')
                if not os.path.exists(rerun_dir): 
                # if already exists, then simulation already re-run- skip this log, move to next
                    os.mkdir(rerun_dir)
                    write_xyz(some_coords, os.path.join(rerun_dir, 'rerun.xyz'))
            else:
                print('No iterations were cycled through!')


    def is_optimisation(self):
        return self.get_runtype() == 'opt'
    
    def is_spec(self):
        return self.get_runtype() == 'sp'

    def is_hessian(self):
        return self.get_runtype() == 'freq'

    def get_data(self):
        """
        Returns the last occurrence of printed energies.
        Currently implemented for DFT only.
        """
        basis = ''
        MP2 = 0.0
        
        for line in self.read():
            if '#P' in line:
                basis = line.split()[1].rsplit('/')[1]
            if re.search('^\sE=\s*-?[0-9]*.[0-9]*', line):
                HF = float(line.split()[1])

        return self.file, self.path, basis, HF, MP2        