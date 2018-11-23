from ..core.results import results

__all__ = ['PsiResults']

class PsiResults(Results):
    """Class defining the results of a PSI4 calculation."""
    def __init__(self, log):
        super().__init__(log)
    
    def __repr__(self):
        return self.log # debugging

    __str__ = __repr__

    def completed(self):
        complete = False
        for line in self.read():
            if 'PSI4 exiting successfully' in line:
                complete = True
        return complete

    def get_runtype(self):
        for line in self.read():
            #regex for energy('mp2') or optimize('scf', dertype='hess') (any number of k-v pairs)
                if re.search("[A-z]*\('[A-z0-9]*'(.?\s*[A-z]*='[A-z]*')*\)", line):
                    if re.search("[A-z]*\('[A-z0-9]*'\)", line): #energy('mp2')
                        return line.split('(')[0]
                    else: #optimize('scf', dertype='hess'......)
                        return line.split('(')[0] #add to this later, using the collect additional data

    def get_hf(self):
        for line in self.read():
            if "Reference Energy" in line:
                res = float(line.split('=')[1].split()[0].strip())
        return res

    def get_e_ss(self):
        for line in self.read():
            if "Same-Spin Energy          =" in line:
                res = float(line.split('=')[1].split()[0].strip())
        return res

    def get_e_os(self):
        for line in self.read():
            if "Opposite-Spin Energy      =" in line:
                res = float(line.split('=')[1].split()[0].strip())
        return res

    def get_mp2(self):
        return self.get_hf() + self.get_e_ss() + self.get_e_os()

    def get_srs(self):
        """Returns SRS-MP2 energy. Note assumes the system contains non-negligible intermolecular
        interactions, such as ionic liquids; the values used in the calculation vary slightly if other systems are used. See J. Chem. Phys. 146, 064108 (2017)"""

        c_os_values = {}
        c_ss_values = {}

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

        return self.get_hf() + c_os * self.get_e_os() + c_ss * self.get_e_ss()
        
    def get_energy(self):
        return self.get_srs()