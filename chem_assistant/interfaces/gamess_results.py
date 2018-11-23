from ..core.utils import write_xyz
from ..core.results import Results

import re

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

    def get_error(self):
        pass

    ################################
    #                              #
    #         IF COMPLETED         #
    #                              #
    ################################


    def get_runtype(self):
        """Returns type of calculation ran"""
        for line in self.read():
            if 'RUNTYP=' in line:
                return line.split()[4].split('=')[1].lower()
        # finds the first instance then breaks out of loop
    
    def get_fmo_level(self):
        """Returns level of FMO calculation ran"""
        for line in self.read():
            if 'NBODY' in line:
              return int(line.split()[-1].split('=')[-1]) # FMO2 or 3
        return 0
        
    def get_equil_coords(self):
        import re
        coords = []
        regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"
        found = False
        for line in self.read(): #generator, so can't use line indexing- returns one line at a time
            if 'EQUILIBRIUM GEOMETRY LOCATED' in line:
                found = True
                # store lines after this if regex found, until INTERNUCLEAR DISTANCES (ANGS.)
            if found:
                if re.search(regex, line):
                    if line.endswith('\n'):
                        coords.append(line[:-1]) # drop newline char
                    else:
                        coords.append(line)
                if 'INTERNUCLEAR DISTANCES (ANGS.)' in line:
                    break 
        write_xyz(coords, 'equil.xyz')   

    ################################
    #                              #
    #      AB INITIO ENERGIES      #
    #                              #
    ################################

    def get_non_fmo(self):
        """Returns energy of the calculation, with the exact definition dependent on the input file parameters. Non-FMO energy"""
        # to store every energy, find line that startwith NSERCH- has iteration and energy on same line
        for line in self.read():
            if 'TOTAL ENERGY =' in line:
                res = float(line.split()[-1])
        return res #last occurrence

    def get_hf(self):
        """Returns Hartree-Fock energy"""
        for line in self.read():
            if 'Euncorr HF(3)=' in line:
                res = float(line.split()[-1])
            elif 'Euncorr HF(2)=' in line:
                res = float(line.split()[-1])
        return res

    def get_mp2(self):
        """Returns MP2 energy"""
        for line in self.read():
            if 'E corr MP2(3)=' in line:
                res = float(line.split()[-1])
            elif 'E corr MP2(2)=' in line:
                res = float(line.split()[-1])
        return res

    def get_srs(self):
        """Returns SRS-MP2 energy. Looks for SCS as this component is set in
the input file using SRS values"""
        for line in self.read():
            if 'E corr SCS(3)=' in line:
                res = float(line.split()[-1])
            elif 'E corr SCS(2)=' in line:
                res = float(line.split()[-1])
        return res

    def get_mp2_corr(self):
        """Returns MP2 correlation energy"""
        for line in self.read():
            if 'Edelta MP2(3)=' in line:
                res = float(line.split()[-1])
            elif 'Edelta MP2(2)=' in line:
                res = float(line.split()[-1])
        return res

    def get_srs_corr(self):
        """Returns SRS-MP2 correlation energy"""
        for line in self.read():
            if 'Edelta SCS(3)=' in line:
                res =  float(line.split()[-1])
            elif 'Edelta SCS(2)=' in line:
                res =  float(line.split()[-1])
        return res

    def get_e_os(self):
        """Calculates the opposite spin component to the interaction energy,
for both MP2 and SRS-MP2 results.

        MP2: Corr = E_OS + E_SS, regardless of basis set.
        SRS: Corr = 1.752*E_OS for cc-pVDZ, 1.64*E_OS for cc-pVTZ

        Other basis sets have coefficients for the same spin parameter, these
need to be accounted for also when the time comes.
        Currently only returns the opposite spin energy found from SRS
calculations due to the greater accuracy when compared to MP2.
        """

        c_os_values['ccd']         = 1.752
        c_os_values['cc-pvdz']     = 1.752
        c_os_values['cct']         = 1.64
        c_os_values['cc-pvtz']     = 1.64
        c_os_values['ccq']         = 1.689
        c_os_values['cc-pvqz']     = 1.689
        c_os_values['accd']        = 1.372
        c_os_values['aug-cc-pvdz'] = 1.372
        c_os_values['acct']        = 1.443
        c_os_values['aug-cc-pvtz'] = 1.443
        c_os_values['accq']        = 1.591
        c_os_values['aug-cc-pvqz'] = 1.591

        c_ss_values['ccd']         = 0
        c_ss_values['cc-pvdz']     = 0
        c_ss_values['cct']         = 0
        c_ss_values['cc-pvtz']     = 0
        c_ss_values['ccq']         = 0
        c_ss_values['cc-pvqz']     = 0
        c_os_values['accd']        = 0
        c_os_values['aug-cc-pvdz'] = 0
        c_os_values['acct']        = 0
        c_os_values['aug-cc-pvtz'] = 0
        c_os_values['accq']        = 0
        c_os_values['aug-cc-pvqz'] = 0

        c_os = c_os_values.get(self.get_basis(), 1.64)
        c_ss = c_ss_values.get(self.get_basis(), 0)

        return self.get_srs_corr() / c_os
        
        # for other basis sets- need the e_ss value- or just solve simultaneously. Come to that when
        # the need arises. This won't work for non-zero c_ss

    def get_e_ss(self):
        if self.get_e_os() is not None:
            return self.get_mp2() - self.get_e_os()
        else:
            return None

    def get_energy(self):
        """Calls different methods depending on the parameters of the input file. The priority of
energies is as follows: FMO3 > FMO2 > Non-FMO.
        The actual call is to either FMO or Non-FMO; the FMO call is then prioritized into FMO3 over
FMO2"""
        level = self.get_fmo_level()
        methods = {0: self.get_non_fmo,
                   2: self.get_srs,
                   3: self.get_srs}

        # need to change in the case of no SRS- need to fall back to get_mp2, and then to
        # get_hf
        return methods.get(level, 2)() # if no FMO, run get_non_fmo()
        

    ################################
    #                              #
    #     VIBRATIONAL ANALYSIS     #
    #                              #
    ################################

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
