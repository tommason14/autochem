import re
import os
from .utils import write_xyz

__all__ = ['Results']

class Results:
    """Base class, only for inheritance"""

    def __init__(self, log):
        self.log = log[2:]
        self.path, self.file = os.path.split(self.log)
        
    def read(self):
        """Memory-efficient reading of large log files, using a generator returning lines as required"""
        with open(self.log, "r") as f:
            for line in f:
                yield line

    def get_error(self):
        print(f'{self.log}: Incomplete calculation')
