from ..core.utils import write_geom_input_for_thermo, write_xyz, eof
from ..core.results import Results

import re
import os
import subprocess
import sys

__all__ = ['GamessResults']

class GamessResults(Results):
    """Class for obtaining results from Gamess simulations. This class requires
    a log file to be read.
    Usage:
        >>> results = Gamess_results(using = 'filename.log')

    Note: assumes that FMO has been used.
    May have to change in future, and write methods for both non-FMO and FMO
    calculations.

    Currently, these methods look for FMO3 data primarily, with an FMO2 fall
back if not found.

    Instances of this class have the following attributes:

    * ``log`` -- filename of the log file of the calculation
    * ``basis`` -- basis set of the calculation, as this class is assumed to be
    * used for ab initio calculations. This attribute may be read from the
    * input file i.e. set as gamess.input.basis = 'CCT')
    * ``coords`` -- coordinates of system, in xyz format
    
    Currently all methods to find energy return the last occurrence of that energy- needs amending to grep every
instance, really. Simple fix; instead of returning values, store in list and return the list, maybe
store the iteration number.
    """
    def __init__(self, log):
        super().__init__(log)

    ################################
    #                              #
    #      CHECK IF COMPLETED      #
    #                              #
    ################################

    def completed(self):
        found = False
        for line in eof(self.log, 0.1):
            if 'EXECUTION OF GAMESS TERMINATED NORMALLY' in line:
                found = True
        return found
        
        ####NEEDS WORK####
        # CURRENTLY IF TERMINATES ABNORMALLY, RESULTS FROM THE CALC
        # ARE NOT RETURNED, EVEN IF THERE

    def memory_error(self):
        print('Memory Error- check allocation before resubmitting')

    def get_error(self):
        super().get_error()
        if self.is_optimisation():
            no_equil = True
            for line in self.read():
                # check for equilibrium coords
                if 'EQUILIBRIUM GEOMETRY LOCATED' in line:
                    no_equil = False    
            if no_equil:
                return 'No equilibrium geometry found- need to resubmit with rerun.xyz'
            else:
                self.memory_error() # last resort

 
    ################################
    #                              #
    #         IF COMPLETED         #
    #                              #
    ################################

    def get_runtype(self):
        """Returns type of calculation ran"""
        for line in self.read():
            if 'RUNTYP=' in line.upper():
                parts = line.split()
                for p in parts:
                    if 'RUNTYP=' in p:
                        return p.split('=')[1].lower()
    
    def get_fmo_level(self):
        """Returns level of FMO calculation ran"""
        for line in self.read():
            if 'NBODY' in line:
              return int(line.split()[-1].split('=')[-1]) # FMO2 or 3
        return 0
        
    def get_equil_coords(self, output = None):
        # find the parent dir for the system, regardless of opt/rerun
        # find the dir with complex/ionic/frags (for frags in subdir) /opt/spec/hess (not frags in
        # subdir)
        # first time that comes up- that's the parent!
        import re
        equil = []
        rerun = []
        regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"
        found_equil = False
        found_some = False
        par_dir = []
        for part in self.path.split('/'):
            if part in ('opt', 'spec', 'hess'):
                break
            else:
                par_dir.append(part)
        MOLECULE_PARENT_DIR = '/'.join(par_dir)
        for line in self.read():
            if 'EQUILIBRIUM GEOMETRY LOCATED' in line:
                found_equil = True
            if 'COORDINATES OF ALL ATOMS ARE (ANGS)' in line: #store every coord list
                found_some = True
                if len(rerun) > 0: # from last run, remove those coords
                    rerun = []
            if found_equil:
                if re.search(regex, line):
                    if line.endswith('\n'):
                        equil.append(line[:-1]) # drop newline char
                    else:
                        equil.append(line)
            if found_some:
                if re.search(regex, line):
                    if line.endswith('\n'):
                        rerun.append(line[:-1]) # drop newline char
                    else:
                        rerun.append(line)
            if line is '\n':
                found_equil = False
                found_some = False
       
        if len(equil) > 0:
            print('found equilibrium!')
            newdir = os.path.join(MOLECULE_PARENT_DIR, 'spec')
            newname = self.basename + '_equil.xyz'
            if not os.path.isdir(newdir):
                os.mkdir(newdir)
            write_xyz(equil, os.path.join(newdir, newname))
        else:
            if len(rerun) > 0:
                print(f'equilibrium not found.\nNeeds resubmitting. Coords stored in {self.path}/rerun/rerun.xyz')
                rerun_dir = os.path.join(self.path, 'rerun')
                if not os.path.exists(rerun_dir): 
                    os.mkdir(rerun_dir)
                write_xyz(rerun, os.path.join(rerun_dir, 'rerun.xyz'))
            else:
                print('No iterations were cycled through!')   
    
    def is_optimisation(self):
        return self.get_runtype() == 'optimize'
    
    def is_spec(self):
        return self.get_runtype() == 'energy'

    def is_hessian(self):
        return self.get_runtype() == 'hessian'

    ################################
    #                              #
    #      AB INITIO ENERGIES      #
    #                              #
    ################################

    
    def fmo_mp2_data(self, mp2_type):
        """
        Returns Hartree Fock and MP2 data.
        Returns the last occurrence of FMO energies 
        (FMO3 given if available, else FMO2), MP2 correlation energies
        and HF energies. Works because FMO3 values are printed 
        after FMO2, and the function returns the last 
        value printed. `mp2_type` should be either 'SCS' or 'MP2', 
        to return the correlated SCS energy, 'E corr SCS', or correlated
        MP2 energies, 'E corr MP2'.
        """
        HF = ''
        MP2 = ''
        for line in eof(self.log, 0.2): # last values only
            if 'Euncorr HF' in line:
                HF = line.split()[-1]
            if f'E corr {mp2_type}' in line:
                MP2 = line.split()[-1]

        HF, MP2 = map(float, (HF, MP2))
        return HF, MP2

    @property
    def total_energy(self):
        """
        Returns last occurrence of total energy.
        """
        total = ''
        for line in self.read():
            if 'TOTAL ENERGY =' in line:
                total = line.split()[-1]
        return float(total)

    @property 
    def basis(self):
        """
        Returns basis set.
        """
        def raw_basis():
            for line in self.read():
                if 'INPUT CARD> $BASIS' in line:
                    return line.split()[-2].split('=')[1]
        
        basis = raw_basis()
        change_basis = {'CCD'  : 'cc-pVDZ',
                        'CCT'  : 'cc-pVTZ',
                        'CCQ'  : 'cc-pVQZ',
                        'aCCD' : 'aug-cc-pVDZ',
                        'aCCT' : 'aug-cc-pVTZ',
                        'aCCQ' : 'aug-cc-pVQZ'}
        return change_basis.get(basis, basis) # if self.basis not in dict, return self.basis
    
    def non_fmo_mp2_data_gas(self):
        """
        Returns value of E(0) as HF, E(2S) as the opposite spin energy and
        E(2T) as same spin energy. Then user can scale energies accordingly.
        """
        HF = ''
        MP2_opp = ''
        MP2_same = ''
        for line in eof(self.log, 0.2):
            line = line.split()
            # print(line)
            if 'E(0)=' in line:
                HF = line[-1]
            if 'E(2S)=' in line:
                MP2_opp = line[-1]
            if 'E(2T)=' in line:
                MP2_same = line[-1]

        HF, MP2_opp, MP2_same = map(float, (HF, MP2_opp, MP2_same))
        return HF, MP2_opp, MP2_same

    def non_fmo_mp2_data_solvent(self):
        """
        Returns value of E(0) as HF, E(MP2) as the MP2 energy. 
        When solvent is added, GAMESS doesn't print the individual correlation
        energy for each spin. Also these energy are for the molecule itself and
        not with the addition of the energy of the solvent. In order to find
        that, search for 'THE P(2) CORRECTED MP2-CPCM ENERGY'.
        """
        HF = ''
        MP2 = ''
        for line in eof(self.log, 0.2):
            line = line.split()
            if 'E(0)=' in line:
                HF = line[-1]
            if 'E(MP2)=' in line:
                MP2 = line[1]

        HF, MP2 = map(float, (HF, MP2))
        return HF, MP2

    @property
    def dft_type(self):
        """
        Returns DFT type (DFTTYP=...)
        """
        for line in self.read():
            if 'DFTTYP' in line:
                line = line.split()
                for val in line:
                    if 'DFTTYP' in val:
                        return val.split('=')[1].upper()

    @property
    def _energy_type(self):
        """
        Returns energy type, i.e. HF, DFT, MP2
        """
        dft = False
        fmo = False
        mp2 = False
        scs = False
        for line in self.read():
            if 'FMO' in line:
                fmo = True
            if 'MPLEVL' in line:
                mp2 = True
            if 'SCS' in line:
                scs = True
            if 'DFT' in line:
                dft = True
            if 'RUN TITLE' in line:
                break # by this point, all data required is specified
        types = {
                 'fmo_scs': (fmo, scs),
                 'fmo_mp2': (fmo, mp2),
                 'fmo_dft': (fmo, dft),
                 'fmo_hf': (fmo,),
                 'scs': (scs,),
                 'mp2': (mp2,),
                 'dft': (dft,)
                }
        for run, possible_runs in types.items():
            if all(r for r in possible_runs):
                return run
        else:
            return 'hf'  

    @property
    def method(self):
        """
        More usable version of self._energy_type.
        i.e. instead of fmo_dft, this method actually returns the 
        dft functional used.
        """
        if 'scs' in self._energy_type:
            return 'Scaled MP2' # scs or srs
        if 'mp2' in self._energy_type and 'scs' not in self._energy_type:
            return 'MP2'
        if 'dft' in self._energy_type:
            return self.dft_type
        if 'hf' in self._energy_type:
            return 'HF'

    @property
    def solvent_calc(self):
        """
        Returns True if the user inputs a $PCM section in the input file.
        No guarantees that the $PCM line will be shown in the copy of the 
        input file at the top though, so instead has to check when the 
        log file reports it.
        """
        for line in self.read():
            if 'INPUT FOR PCM SOLVATION CALCULATION' in line:
                return True
            if 'BEGINNING GEOMETRY SEARCH' in line:
                break
        return False

    def get_data(self):
        """
        Returns HF, DFT, MP2 (non-FMO and FMO) values.
        """
        if self._energy_type == 'fmo_scs':
            HF, MP2 = self.fmo_mp2_data('SCS')
            MP2_opp = 'NA'
            MP2_same = 'NA'
        elif self._energy_type == 'fmo_mp2':
            HF, MP2 = self.fmo_mp2_data('MP2')
            MP2_opp = 'NA'
            MP2_same = 'NA'
        elif self._energy_type in ('mp2', 'scs'):
            if not self.solvent_calc:
                HF, MP2_opp, MP2_same = self.non_fmo_mp2_data_gas()
                MP2 = 'NA'
            else:
                HF, MP2 = self.non_fmo_mp2_data_solvent()
                MP2_opp = 'NA'
                MP2_same = 'NA'
        elif 'dft' in self._energy_type or 'hf' in self._energy_type:
            HF = self.total_energy
            MP2 = 'NA'
            MP2_opp = 'NA'
            MP2_same = 'NA'
            
        return self.file, self.path, self.method, self.basis, HF, MP2, MP2_opp, MP2_same    

    @property
    def multiplicity(self):
        for line in self.read():
            if 'SPIN MULTIPLICITY' in line.upper(): # sometimes prints lower case
                return int(line.split()[-1])

    #################
    ### HOMO-LUMO ###
    #################

    @property
    def num_orbitals_occupied(self):
        for line in self.read():
            if 'ORBITALS ARE OCCUPIED' in line:
                return int(line.split()[0])

    def _homo_lumo(self):
        """
        Finds HOMO/LUMO orbitals.
        """
        found = False
        orbital_energies = []
        for line in self.read():
            if 'EIGENVECTORS' in line:
                found = True
            if 'CPU' in line:
                found = False
            if found:
                if re.search('^(\s+-?[0-9]+.[0-9]+){1,5}$', line):
                    orbital_energies += [float(i) for i in line.split()]
                    # save time
                    if len(orbital_energies) > self.num_orbitals_occupied + 1:
                        break
        homo = orbital_energies[self.num_orbitals_occupied - 1]
        lumo = orbital_energies[self.num_orbitals_occupied]
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

            
    ################################
    #                              #
    #     VIBRATIONAL ANALYSIS     #
    #                              #
    ################################

    # create and visualise freq modes/ IR spectra.
    # If I ever fit multiple lorentzians to the peaks, won't need molden
    # anymore!
    # Then this will become useful and could replace the gamess_to_molden script.

    def vib_get_geom(self):
        pass
    
    @property
    def frequencies(self):
        """
        Returns vibrations, with translations/rotations included.
        Checks output below this line:
        'MODE FREQ(CM**-1)  SYMMETRY  RED. MASS  IR INTENS.'
        """
        vibs = []
        regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
        for line in self.eof(0.2): # last 20 % of file
            if re.search(regex, line):
                vibs.append(float(line.split()[1]))
        return vibs

    @property
    def intensities(self):
        """
        Returns intensities of all vibrations
        Checks output below this line:
        'MODE FREQ(CM**-1)  SYMMETRY  RED. MASS  IR INTENS.'
        """ 
        ints = []
        regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
        for line in self.eof(0.2): # last 20 % of file
            if re.search(regex, line):
                ints.append(float(line.split()[-1]))
        return ints

    def write_initial_geom_for_thermo(self):
        """Parses GAMESS inputs for the initial geometry"""
        atoms = []
        regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"
        inp = file[:-3] + 'inp'
        if inp not in os.listdir('.'):
            sys.exit(f'Need an input file in the same directory as {self.log}')
        for line in read_file(inp):
            if re.search(regex, line):
                sym, _, x, y, z = line.split()
                x, y, z = map(float, (x, y, z))
                atoms.append(Atom(symbol = sym, coords = (x, y, z)))
        write_geom_input_for_thermo(atoms)

