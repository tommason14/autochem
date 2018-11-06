#!/usr/bin/env python3

from molecule import Molecule

mol = Molecule('../test_xyz/c1mim-nh3.xyz')
mol.separate()
mol.gamess_format()