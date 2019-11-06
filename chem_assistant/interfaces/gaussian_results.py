from ..core.utils import read_file, write_xyz
from ..core.results import Results
from ..core.periodic_table import PeriodicTable as PT
from ..core.atom import Atom

import re
import os
import subprocess

__all__ = ['GaussianResults']

class GaussianResults(Results):
    """
    Class for obtaining results from Gaussian simulations. This class requires a log file to be read.
    """

    def __init__(self, log):
        super().__init__(log)

    ################################
    #                              #
    #      CHECK IF COMPLETED      #
    #                              #
    ################################

    def get_runtype(self):
        """
        Can have multiple runs in one file i.e. opt freq
        Looks for one or both of 'opt' and 'freq'.
        If not found, defaults to single point calculation.
        """
        opt = False
        freq = False
        for p in self.user_commands.split():
            if 'opt' in p:
                opt = True
            if 'freq' in p:
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
        found = False
        for line in self.eof(0.01):
            if 'Normal termination' in line:
                return True
        return False

    def _calcall_equil_coords(self):
        """
        For opt=(...calcall), coords printed differently
        """
        regex = "^(\s+[0-9]+){3}(\s+-?[0-9]{1,3}.[0-9]+){3}$"
        coords = []
        found_some = False
        found_equil = False
        for line in self.read():
            if 'Optimization completed' in line:
                found_equil = True
            if 'Standard orientation' in line:
                found_some = True
                if len(coords) > 0:
                    coords = [] # after each iter
            if 'Rotational constants' in line:
                found_some = False
            # check that opt has finished (if opt freq ran)
            if re.search('Freq$', line):
                break
            if found_some:
                if re.search(regex, line):
                    _, atnum, _, x, y, z = line.split()
                    atnum, x, y, z = map(float, (atnum, x, y, z))
                    coords.append(Atom(atnum=atnum, coords = (x, y, z)))
        if found_equil:
            return coords, 'equil'
        else:
            return coords, 'rerun'

    def get_equil_coords(self, output = None):
        """
        Returns either the equilibrium coordinates or coordinates for the rerun.
        """
        # with calcall, things are different, collect standard orientation
        # and if optimization completed found, they are equil coords
        # else rerun
        calcall = True if 'calcall' else False in self.user_commands
        found_equil = False
        found_some = False
        regex = "^(\s+[0-9]+){3}(\s+-?[0-9]{1,3}.[0-9]+){3}$"
        coords = []
        some_coords = []
        par_dir = []
        for part in self.path.split('/'):
            if part not in ('opt', 'spec', 'hess'):
                par_dir.append(part)
            else:
                break
        MOLECULE_PARENT_DIR = '/'.join(par_dir)
        if not calcall:
            for line in self.read():
                if 'Optimization completed' in line:
                    found_equil = True
                if 'Standard orientation' in line:
                    found_some = True
                    if len(some_coords) > 0: # from last run, remove those coords
                        some_coords = []
                if found_equil:
                    if re.search(regex, line):
                        _, atnum, _, x, y, z = line.split()
                        atnum, x, y, z = map(float, (atnum, x, y, z))
                        coords.append(Atom(atnum=atnum, coords = (x, y, z)))
                if found_some:
                    if re.search(regex, line):
                        _, atnum, _, x, y, z = line.split()
                        atnum, x, y, z = map(float, (atnum, x, y, z))
                        some_coords.append(Atom(atnum=atnum, coords = (x, y, z)))
                if 'Rotational constants' in line:
                    found_some = False
                if 'Distance matrix (angstroms)' in line:
                    found_equil = False
        
        if calcall:
            coords_found, coord_type = self._calcall_equil_coords()
            if coord_type == 'equil':
                coords = coords_found
            else:
                some_coords = coords_found
        if len(coords) > 0:
            print('found!')
            write_xyz(coords, os.path.join(MOLECULE_PARENT_DIR, f'{self.title}_equil.xyz'))
        else:
            if len(some_coords) > 0:
                print(f'not found.\nNeeds resubmitting. Coords stored in {self.path}/rerun.xyz')
                write_xyz(some_coords, os.path.join(MOLECULE_PARENT_DIR, 'rerun.xyz'))
            else:
                print('No iterations were cycled through!')
        
    @property
    def title(self):
        """
        Finds title of the job. Printed like so:
        ...
        -----
        title
        -----
        Symbolic Z-matrix
        ...
        """
        lines = []
        found = False
        for line in self.read():
            if 'Leave Link' in line:
                found = True
                continue 
            if 'Symbolic' in line:
                break
            if found:
                lines.append(line)
        return lines[-2].strip()

    @property
    def user_commands(self):
        """
        Returns the #P line of the input file.
        """
        for line in self.read():
            if re.search('^\s*?#P', line):
                return line.lower()

    @property
    def method(self):
        """
        Returns energy type. For example, for HF/cc-pVTZ, returns hf. 
        For wB97xD/aug-cc-pVDZ, returns wb97xd.
        """
        return self.user_commands.split('/')[0].split()[-1].upper()

    @property
    def basis(self):
        """
        Returns basis set. For example, for HF/cc-pVTZ, returns cc-pvtz. 
        For wB97xD/aug-cc-pVDZ, returns aug-cc-pVDZ.
        """
        basis = self.user_commands.split('/')[1].split()[0]
        # turn cc-pvtz into cc-pVTZ
        if 'cc' in basis:
            *rest, last = basis.split('-')
            last = last[:-3] + last[-3:].upper()
            basis = '-'.join(rest + [last])
        return basis

    def is_optimisation(self):
        return 'opt' in self.get_runtype()

    def is_spec(self):
        return self.get_runtype() == 'spec'

    def is_hessian(self):
        return 'freq' in self.get_runtype()

    @property
    def hf_energy(self):
        """
        Returns last occurrence of Hartree-Fock energy.
        """
        HF = ''
        for line in self.read():
            if re.search('^\sE=\s*-?[0-9]*.[0-9]*', line):
                HF = line.split()[1]
        return float(HF)
    
    @property
    def mp2_energy(self):
        """
        Returns last occurrence of MP2 energy.
        """
        # needs filling
        return 0.0

    @property
    def dft_energy(self):
        """
        Returns last occurrence of DFT energy.
        """
        dft = ''
        for line in self.read():
            if 'SCF Done' in line:
                dft = line.split()[4]
        return float(dft)

    def get_data(self):
        """
        Returns the last occurrence of printed energies. Negate energy types, and if the energy type
        is not found, assumed to be DFT. 
        Must return file, path, method, basis, hf/dft, mp2/srs, mp2_opp, mp2_same
        in the order, so 'NA' values are there to satisfy that criteria.
        """
        if self.method == 'hf':
            return self.file, self.path, self.method, self.basis, self.hf_energy, 'NA', 'NA', 'NA'
        elif self.method == 'mp2':
            return self.file, self.path, self.method, self.basis, self.hf_energy, self.mp2_energy, 'NA', 'NA'
        else:
            return self.file, self.path, self.method, self.basis, self.dft_energy, 'NA', 'NA', 'NA'

    @property
    def multiplicity(self):
        """
        Returns multiplicity from the symbolic z-matrix section.
        """
        for line in self.read():
            if 'Charge' in line and 'Multiplicity' in line:
                return int(line.split()[-1])
    
    def _homo_lumo(self):
        """
        Finds HOMO/LUMO orbitals.
        """
        occupied = []
        lumo = ''
        for line in self.read():
            if 'Alpha  occ. eigenvalues' in line:
                occupied += line.split()[4:]
            if 'Alpha virt. eigenvalues' in line:
                lumo = line.split()[4] 
                break
        homo = occupied[-1]
        homo, lumo = map(float, (homo, lumo))
        return homo, lumo

    def _homo_lumo_gap(self):
        hartrees_to_eV = 27.21
        homo, lumo = self._homo_lumo()
        homo, lumo = map(lambda x: x * hartrees_to_eV, (homo, lumo))
        gap = lumo - homo
        return homo, lumo, gap

    @property
    def homo_lumo_info(self):
        """
        Prints the HOMO-LUMO gap. Finds SOMO-LUMO if multiplicity is 2.
        Returns `self.multiplicity`, SOMO/HOMO (eV), LUMO (eV) and the gap (eV).
        Note that ⍺ orbitals are selected, and the ⍺ and β electrons give
        different homo-lumo gaps.
        """

        if self.multiplicity == 1:
            homo, lumo, gap = self._homo_lumo_gap()
            transition = 'HOMO-LUMO'
        elif self.multiplicity == 2:
            homo, lumo, gap = self._homo_lumo_gap() # here homo is somo
            transition = 'SOMO-LUMO'
        else:
            print(f'Error: Only singlet/doublet multiplicities have been accounted for. Ignoring {self.log}')
            
        return {'File': self.file,
                'Path': self.path,
                'Multiplicity': self.multiplicity, 
                'Transition': transition, 
                'HOMO/SOMO (eV)': homo, 
                'LUMO (eV)': lumo, 
                'Gap (eV)': gap}

