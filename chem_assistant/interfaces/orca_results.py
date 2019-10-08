from ..core.utils import read_file, write_xyz
from ..core.results import Results
from ..core.periodic_table import PeriodicTable as PT
from ..core.atom import Atom

import re
import os
import subprocess

__all__ = ['OrcaResults']

class OrcaResults(Results):
    """
    Class for obtaining results from Orca simulations. This class     
    requires a log file to be read.
    """

    def __init__(self, log):
        super().__init__(log)

    def __repr__(self):
        return self.log

    __str__ = __repr__

    def get_runtype(self):
        """
        Can have multiple runs in one file i.e. opt freq
        Looks for one or both of 'opt' and 'freq'.
        If not found, defaults to single point calculation.
        """
        opt = False
        freq = False
        if 'opt' in self.user_commands:
            opt = True
        if 'freq' in self.user_commands:
            freq = True
        if opt and not freq:
            return 'opt'
        if opt and freq: 
            return 'opt-freq'
        if freq and not opt:
            return 'freq'
        if not opt and not freq:
            return 'spec'

    def completed(self):
        for line in self.eof(0.05):
            if '****ORCA TERMINATED NORMALLY****' in line:
                return True
        return False

    # def get_equil_coords(self, output = None):
    #     found_equil = False
    #     found_some = False
    #     regex = "^(\s*[0-9]*){3}(\s*-?[0-9]{1,3}.[0-9]*){3}$"
    #     coords = []
    #     some_coords = []
    #     par_dir = []
    #     # print(self.log)
    #     for part in self.path.split('/'):
    #         if part not in ('opt', 'spec', 'hess'):
    #             par_dir.append(part)
    #         else:
    #             break
    #     MOLECULE_PARENT_DIR = '/'.join(par_dir)
    #     for line in self.read():
    #         if 'Optimization completed' in line:
    #             found_equil = True
    #         if 'Standard orientation' in line:
    #             found_some = True
    #             if len(some_coords) > 0: # from last run, remove those coords
    #                 some_coords = []
    #         if found_equil:
    #             if re.search(regex, line):
    #                 _, atnum, _, x, y, z = line.split()
    #                 atnum, x, y, z = map(float, (atnum, x, y, z))
    #                 sym = PT.get_symbol(atnum)
    #                 coords.append(Atom(sym, coords = (x, y, z)))
    #         if found_some:
    #             if re.search(regex, line):
    #                 _, atnum, _, x, y, z = line.split()
    #                 atnum, x, y, z = map(float, (atnum, x, y, z))
    #                 sym = PT.get_symbol(atnum)
    #                 some_coords.append(Atom(sym, coords = (x, y, z)))
    #         if 'Rotational constants' in line:
    #             found_some = False
    #         if 'Distance matrix (angstroms)' in line:
    #             found_equil = False
    #     if len(coords) > 0:
    #         print('found!')
    #         write_xyz(coords, os.path.join(MOLECULE_PARENT_DIR, 'equil.xyz'))
    #     else:
    #         if len(some_coords) > 0:
    #             print(f'not found.\nNeeds resubmitting. Coords stored in {self.path}/rerun.xyz')
    #             write_xyz(some_coords, os.path.join(MOLECULE_PARENT_DIR, 'rerun.xyz'))
    #             # rerun_dir = os.path.join(self.path, 'rerun')
    #             # if not os.path.exists(rerun_dir):
    #             # # if already exists, then simulation already re-run- skip this log, move to next
    #             #     os.mkdir(rerun_dir)
    #             #     write_xyz(some_coords, os.path.join(rerun_dir, 'rerun.xyz'))
    #         else:
    #             print('No iterations were cycled through!')


    def is_optimisation(self):
        return 'opt' in self.get_runtype()

    def is_spec(self):
        return self.get_runtype() == 'spec'

    def is_hessian(self):
        return 'freq' in self.get_runtype()
    
    @property
    def user_commands(self):
        """
        Returns the ! ... line of the input file.
        """
        for line in self.read():
            if '> !' in line:
                return line.lower()

    @property
    def energy_type(self):
        """
        Returns energy type.
        Currently parses the line of ! .... commands for specific values.
        """
        types = ('hf', 'mp2', 'cc')

        for t in types:
            if t in self.user_commands:
                return t
        # return t for t in types if t in self.user_commands

    @property
    def basis(self):
        """
        Returns basis set.
        """
        for line in self.read():
            if 'Your calculation utilizes the basis:' in line:
                return line.split()[-1]

    @property
    def total_energy(self):
        """
        Returns total energy, printed for scf calculations.
        """
        for line in self.read():
            if 'Total Energy       :' in line:
                return float(line.split()[3].strip())
 
    def scf_data(self):
        """
        Return data for scf calculations.
        Note the NAs returned are because of no MP2 spin parameters.
        """
        return self.file, self.path, self.basis, self.total_energy, 'NA', 'NA'


    def get_data(self):
        """
        Returns the last occurrence of printed energies.
        Currently implemented for DFT only.
        """

    def get_data(self):
        """
        Returns job data: filename, filepath, basis set, HF/DFT energy, and MP2 opposite
        and same spin parameters if relevant.
        """
        if self.energy_type == 'hf':
            return self.scf_data()

        # elif self.energy_type == 'mp2':
            # return self.mp2_data()
    
