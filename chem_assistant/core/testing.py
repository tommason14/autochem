#!/usr/bin/env python3

from gamess import GamessJob
# from molecule import Molecule
from settings import Settings


s = Settings()
s.input.contrl.runtyp = 'hessian'
j = GamessJob(using = '../test_xyz/ch_ac.xyz', settings = s, fmo=True, frags_in_subdir=True) #automatically creates input


