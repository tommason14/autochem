import re
import os
from .utils import write_xyz, eof

__all__ = ['Results']

class Results:
    """Base class, only for inheritance"""

    def __init__(self, log):
        # self.log = log[2:]
        self.log = log
        self.path, self.file = os.path.split(self.log)
        self.basename = self.file.split('.')[0]
        self.abspath = os.path.abspath(log)
        self.parent_dir = self.abspath.split('/')[-2]
        
    def read(self):
        """
        Memory-efficient reading of large log files, using a generator 
        returning lines as required
        """
        with open(self.log, "r") as f:
            for line in f:
                yield line

    def get_error(self):
        print(f'{self.log}: Incomplete calculation')

    def eof(self, percentage):
        """
        Include percentage as decimal.
        i.e. self.eof(0.05) returns the last 5% of the file
        """
        return eof(self.abspath, percentage) 
