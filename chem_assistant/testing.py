#!/usr/bin/env python3

# from molecule import Molecule

# mol = Molecule('../test_xyz/ch_ac.xyz')
# mol.separate()
# mol.gamess_format()
from input import Gamess_input

gamess = Gamess_input(using = '../test_xyz/water.xyz', fmo = True)