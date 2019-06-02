from ..core.results import Results

import re

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
            if 'Psi4 exiting successfully' in line:
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

    def is_optimisation(self):
        return self.get_runtype() == 'optimize'
    
    def is_spec(self):
        return self.get_runtype() == 'energy'

    def is_hessian(self):
        return self.get_runtype() == 'frequency'

    def get_data(self):
        HF = ''
        opp = ''
        same = ''
        basis = ''
        MP2 = '' 
        # with open(filepath, "r") as f:
        #     for line in f.readlines():
        for line in self.read():
            if re.search('basis\s\w*(\-?\w*){1,2}$', line):
                basis = line.split()[-1]
            if 'Reference Energy          =' in line:
                HF = float(line.split('=')[1].split()[0].strip())
            elif 'Same-Spin Energy          =' in line:
                same = float(line.split('=')[1].split()[0].strip())
            elif 'Opposite-Spin Energy      =' in line:
                opp = float(line.split('=')[1].split()[0].strip())
        
        SRS = {}
        SRS['c_os'] = {}
        SRS['c_os']['ccd']         = 1.752
        SRS['c_os']['cc-pvdz']     = 1.752
        SRS['c_os']['cct']         = 1.64
        SRS['c_os']['cc-pvtz']     = 1.64
        SRS['c_os']['ccq']         = 1.689
        SRS['c_os']['cc-pvqz']     = 1.689
        SRS['c_os']['accd']        = 1.372
        SRS['c_os']['aug-cc-pvdz'] = 1.372
        SRS['c_os']['acct']        = 1.443
        SRS['c_os']['aug-cc-pvtz'] = 1.443
        SRS['c_os']['accq']        = 1.591
        SRS['c_os']['aug-cc-pvqz'] = 1.591

        SRS['c_ss'] = {}
        SRS['c_ss']['ccd']         = 0
        SRS['c_ss']['cc-pvdz']     = 0
        SRS['c_ss']['cct']         = 0
        SRS['c_ss']['cc-pvtz']     = 0
        SRS['c_ss']['ccq']         = 0
        SRS['c_ss']['cc-pvqz']     = 0
        SRS['c_ss']['accd']        = 0
        SRS['c_ss']['aug-cc-pvdz'] = 0
        SRS['c_ss']['acct']        = 0
        SRS['c_ss']['aug-cc-pvtz'] = 0
        SRS['c_ss']['accq']        = 0
        SRS['c_ss']['aug-cc-pvqz'] = 0


        c_os = SRS['c_os'][basis.lower()]
        c_ss = SRS['c_ss'][basis.lower()]
        
        MP2 = HF + c_os * opp + c_ss * same
        
        return self.file, self.path, basis, HF, MP2

