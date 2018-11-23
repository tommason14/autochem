import re
from .utils import write_xyz

__all__ = ['Results']

class Results:
    """Base class, only for inheritance"""

    def __init__(self, log):
        self.log = log
        # if coords is not None:
        #     self.input_coords = input_coords
        #     # coords needed for hessian calcs

    def read(self):
        """Memory-efficient reading of large log files, using a generator returning lines as required"""
        with open(self.log, "r") as f:
            for line in f.readlines():
                yield line

    def get_type(self):
        for line in self.read():
            if 'PSI4' in line:
                return 'psi4'
            elif 'GAMESS' in line:
                return 'gamess'

    def get_basis(self):
        for line in self.read():
            # include file type as a precaution ('basis' could appear in a gamess file)
            if self.get_type() == 'psi4':
                if 'basis' in line.lower():
                    return line.split()[-1].lower()
            elif self.get_type() == 'gamess':
                if 'gbasis' in line.lower(): #gamess data can be written as upper or lowercase
                    return line.split()[-2].split('=')[-1].lower()