import os
import subprocess

__all__ = ['Supercomp']

class Supercomp:
    """Detects the supercomputer in use when called. In reality, using a class
    may be overkill here. However, this allows for extensibility when required.
    A instance of this class can be concatenated with strings- used internally 
    when writing job files:
        >>> sc = Supercomp()
        >>> print('gamess_' + sc) # gamess_rjn or gamess_mgs   
    """

    def __init__(self):
        cases = {'raijin': 'rjn',
                 'magnus': 'mgs',
                 'nfs': 'gaia',
                 'm3': 'mas',
                 'monarch': 'mon',
                 'stampede': 'stm'
                 }
        hostname = subprocess.run('hostname', encoding = 'utf-8',
                   stdout = subprocess.PIPE).stdout[:-1]
        for key in cases:
            if key in hostname:
                self.sc = cases[key]
                break
            else:
                self.sc = 'stm' 

    def __repr__(self):
        return str(self.sc)

    def __add__(self, other):
        """String concatenation"""
        if isinstance(other, str):
            return self.sc + other
        else:
            raise TypeError('Can only concatenate Supercomp instances with strings. No addition to floats/integers')
    
    def __radd__(self, other):
        """String concatenation"""
        if isinstance(other, str):
            return other + self.sc
        else:
            raise TypeError('Can only concatenate Supercomp instances with strings. No addition to floats/integers')
   
    @property
    def supercomp(self):
        return self.sc
    
    __str__ = __repr__
