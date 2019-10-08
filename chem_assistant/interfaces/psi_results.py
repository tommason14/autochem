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
            if 'exiting successfully' in line:
                complete = True
        return complete

    def get_runtype(self):
        """
        Returns runtype. For example, for MP2 single points, the line `energy('mp2')` is used. 
        This method returns the string 'energy'.
        """
        for line in self.read():
            # need regex for energy('mp2') or optimize('scf', dertype='hess') (any number of k-v pairs)
                if re.search("[A-z]*\('[A-z0-9]*'(.?\s*[A-z]*='[A-z]*')*\)", line):
                    if re.search("[A-z]*\('[A-z0-9]*'\)", line): #energy('mp2')
                        return line.split('(')[0]
                    else: #optimize('scf', dertype='hess'......)
                        return line.split('(')[0] #add to this later, using the collect additional data

    @property
    def energy_type(self):
        """
        Returns energy type. For example, for MP2 single points, the line `energy('mp2')` is used. 
        This method returns the string 'mp2'.
        """
        for line in self.read():
            # need regex for energy('mp2') or optimize('scf', dertype='hess') (any number of k-v pairs)
                if re.search("[A-z]*\('[A-z0-9]*'(.?\s*[A-z]*='[A-z]*')*\)", line):
                    if re.search("[A-z]*\('[A-z0-9]*'\)", line): #energy('mp2')
                        return re.search("[A-z]*\('([A-z0-9]*)'\)", line).group(1)
                    # else: #optimize('scf', dertype='hess'......)
                        # return line.split('(')[0] #add to this later, 

    def is_optimisation(self):
        return self.get_runtype() == 'optimize'
    
    def is_spec(self):
        return self.get_runtype() == 'energy'

    def is_hessian(self):
        return self.get_runtype() == 'frequency'

    def get_data(self):
        """
        Returns job data: filename, filepath, basis set, HF/DFT energy, and MP2 opposite
        and same spin parameters if relevant.
        """  
        if self.energy_type == 'scf':
            return self.scf_data()
        
        elif self.energy_type == 'mp2':
            return self.mp2_data() 

    @property
    def basis(self):
        """
        Returns basis set.
        """ 
        for line in self.read():
            if re.search('basis\s\w*(\-?\w*){1,2}$', line):
                return line.split()[-1]

    @property
    def total_energy(self):
        """
        Returns total energy, printed for scf calculations.
        """
        total=''
        for line in self.read():
            if 'Total Energy =' in line:
                total = float(line.split('=')[1].strip())
        return total

    def scf_data(self):
        """
        Return data for scf calculations. 
        Note the NAs returned are because of no MP2 data.
        """
        return self.file, self.path, self.basis, self.total_energy, 'NA', 'NA', 'NA'

    @property
    def hf_energy_for_mp2(self):
        """
        Returns 'reference energy' from MP2 calculations.
        """
        HF=''
        for line in self.read():
            if 'Reference Energy          =' in line:
                HF = float(line.split('=')[1].split()[0].strip())
        return HF

    @property
    def mp2_opp(self):
        """
        Returns MP2 opposite spin energy.
        """
        opp=''
        for line in self.read():
            if 'Opposite-Spin Energy      =' in line:
                opp = float(line.split('=')[1].split()[0].strip())
        return opp 
        
    @property
    def mp2_same(self):
        """
        Returns MP2 same spin energy.
        """
        same=''
        for line in self.read():
            if 'Same-Spin Energy          =' in line:
                same = float(line.split('=')[1].split()[0].strip())
        return same
    
    def mp2_data(self):
        """
        Returns data for MP2 calculations: filename, filepath, 
        basis set, hf energy, opp spin energy, same spin energy.
        No MP2 data is returned, but is calculated instead from the 
        HF and MP2 correlation energies.
        """
        
        return (self.file, self.path, self.basis, 
        self.hf_energy_for_mp2, 'NA', self.mp2_opp, self.mp2_same)
