import re
import os
from .utils import write_xyz

__all__ = ['Results']

class Results:
    """Base class, only for inheritance"""

    def __init__(self, log):
        self.log = log
        self.filepath = os.path.split(self.log)[0]
        self.abspath = os.path.split(os.path.abspath(self.log))[0]
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
        # include file type as a precaution ('basis' could appear in a gamess file)
        if self.get_type() == 'psi4':
            for line in self.read():
                if 'basis' in line.lower():
                    return line.split()[-1].lower()
        elif self.get_type() == 'gamess':
            for line in self.read():
                if 'gbasis' in line.lower() and 'input card' in line.lower(): #gamess data can be written as upper or lowercase
                    return line.split()[-2].split('=')[-1].lower()
        else:
            return None

    def get_error(self):
        print(f'{self.log}: Incomplete calculation')
