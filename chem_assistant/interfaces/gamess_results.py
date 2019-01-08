from ..core.utils import write_xyz
from ..core.results import Results

import re
import os
import subprocess

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

    def __repr__(self):
        return self.log # debugging

    __str__ = __repr__

    def completed(self):
        found = False
        for line in self.read():
            if 'EXECUTION OF GAMESS TERMINATED NORMALLY' in line:
                found = True
        return found

    # call like this:
    # if not self.completed():
    #     self.get_error() 
    # cutoff in middle of run with no explanation- memory error
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
        # finds the first instance then breaks out of loop
    
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
        print(self.log)
        for part in self.path.split('/'):
            if part not in ('opt', 'spec', 'hess'):
                par_dir.append(part)
            else:
                break
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
            print('found!')
            # find parent dir- if multiple reruns, then don't know how far up tree to go to find
            # parent
            write_xyz(equil, os.path.join(MOLECULE_PARENT_DIR, 'equil.xyz'))
            # self.create_spec_after_opt()
        else:
            if len(rerun) > 0:
                print(f'not found.\nNeeds resubmitting. Coords stored in {self.path}/rerun/rerun.xyz')
                # make new dir, and copy over input, replacing any coords
                # TODO: Refactor this, but cba right now
                rerun_dir = os.path.join(self.path, 'rerun')
                if not os.path.exists(rerun_dir): 
                # if already exists, then simulation already re-run- skip this log, move to next
                    os.mkdir(rerun_dir)
                    write_xyz(rerun, os.path.join(rerun_dir, 'rerun.xyz'))
                    basename, ext = self.file.split('.')
                    inp = basename + '.inp'
                    job = basename + '.job'
                    orig_inp = os.path.join(self.path, inp) # path of the log file
                    orig_job = os.path.join(self.path, job)
                    rerun_inp = os.path.join(rerun_dir, 'rerun.inp')
                    rerun_job = os.path.join(rerun_dir, 'rerun.job')
                    print(orig_inp)                    
                    print(orig_job) 
                    print(rerun_inp)
                    print(rerun_job)

                    os.system(f'sed "s/{basename}/rerun/g" {orig_job} >> {rerun_job}') # opt.inp --> rerun.inp
                    # os.system(f'sed "s/{self.file}/rerun.{ext}/g" rerun/rerun_job') # opt.log --> rerun.log
                    # parse original inp and add new coords
                    rerun_inp_file = []
                    with open(orig_inp, "r") as f:
                        for line in f.readlines():
                            if re.search(regex, line):
                                break
                            else:
                                rerun_inp_file.append(line)
                    for line in rerun: #add coords
                        rerun_inp_file.append(line + '\n') 
                    rerun_inp_file.append(' $END')
                    with open(rerun_inp, "w") as f:
                        for line in rerun_inp_file:
                            f.write(line)
            else:
                print('no iterations were cycled through!')   
    
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

    def calc_type(self):
        fmo = False
        mp2 = False
        scs = False
        # fmo, mp2/srs, hf, dft?
        # with open(filepath, "r") as f:
        #     for line in f.readlines():
        for line in self.read():
            if 'FMO' in line:
                fmo = True
            elif 'MPLEVL' in line:
                mp2 = True
            elif 'SCS' in line:
                scs = True
            elif 'RUN TITLE' in line:
                break # by this point, all data required is specified
        return fmo, mp2, scs

    def mp2_data(self, mp2_type):
        basis = ''
        HF = ''
        MP2 = ''
        # with open(filepath, "r") as f:
        #     for line in f.readlines():
        for line in self.read():
            if 'Euncorr HF' in line:
                HF = float(line.split()[-1])
            if 'INPUT CARD> $BASIS' in line:
                basis = line.split()[-2].split('=')[1]
            if f'E corr {mp2_type}' in line:
                MP2 = float(line.split()[-1])
        return basis, HF, MP2

    def get_data(self):
        """Returns the last occurrence of FMO energies (FMO3 given if available, else FMO2), SRS and HF energies"""
        fmo, mp2, scs = self.calc_type()
        if fmo and scs:
            basis, HF, MP2 = self.mp2_data('SCS')
        elif fmo and mp2 and not scs:
            basis, HF, MP2 = self.mp2_data('MP2')
        elif not fmo:
            val = ''
            basis = ''
            HF = ''
            MP2 = ''
            
            # with open(filepath, "r") as f:
            #     for line in f.readlines():
            for line in self.read():
                if 'INPUT CARD> $BASIS' in line:
                    basis = line.split()[-2].split('=')[1]
                if 'TOTAL ENERGY =' in line:
                    val = float(line.split()[-1])
            # What is total energy? SRS/MP2 or HF/DFT or something else?
            if scs or mp2:
                MP2 = val
            else:
                HF = val

        # more readable basis set
        change_basis = {'CCD'  : 'cc-pVDZ',
                        'CCT'  : 'cc-pVTZ',
                        'CCQ'  : 'cc-pVQZ',
                        'aCCD' : 'aug-cc-pVDZ',
                        'aCCT' : 'aug-cc-pVTZ',
                        'aCCQ' : 'aug-cc-pVQZ'}
        basis = change_basis.get(basis, basis) # default is the current value

        return self.file, self.path, basis, HF, MP2        


            
    ################################
    #                              #
    #     VIBRATIONAL ANALYSIS     #
    #                              #
    ################################

    # create and visualise freq modes/ IR spectra    

    def vib_get_geom(self):
        pass
    
    def vib_get_vibs(self):
        vibs = {}
        pass
        # vibs[mode] = [list of vibs for each atom]
        # one set for each mode
    
    def vib_get_intensity(self):
        ints = {}
        pass
        # ints[mode] = [list of ints for each atom]

    def vib_get_coords(self):
        coords = {}
        pass
        # coords[mode] = [list of coords, one coord for each atom]
