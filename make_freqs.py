#!/usr/bin/env python3
from chem_assistant import GamessJob, Settings
import os

par = os.getcwd()
s = Settings()
s.input.contrl.runtyp = 'hessian'
s.input.force.method = 'seminum'
for path, _, files in os.walk('.'):
    for file in files:
        if file == 'equil.xyz' and 'cm5' not in path:
            os.chdir(path)
            print(os.getcwd())
            GamessJob(fmo = True, frags_in_subdir = True, using = 'equil.xyz', settings = s,
is_complex = True)
            os.chdir(par)

