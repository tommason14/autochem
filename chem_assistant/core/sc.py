import os

__all__ = ['Supercomp']

class Supercomp:
    """Detects the supercomputer in use when called. In reality, using a class may be overkill here. However, this allows for extensibility when required.
    A instance of this class can be concatenated with strings- used internally when writing job files:
        >>> sc = Supercomp()
        >>> print('gamess_' + sc) # gamess_rjn or gamess_mgs
"
    Note: On Magnus, the script assumes it is used in the 'scratch' directory- 'scratch' must appear somewhere in the current path. If not, this class will produce the wrong result by assuming use on a local machine and not on a remote file system.
    """

    def __init__(self):
        cwd = os.getcwd()
        cases = {'565': 'rjn',
                 'scratch': 'mgs',
                 'nfs': 'gaia',
                 }
        for key in cases:
            if key in cwd:
                self.sc = cases[key]
                break
            else:
                self.sc = 'rjn' #default to raijin

    def __repr__(self):
        return f"{self.sc}"

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
    
    __str__ = __repr__
