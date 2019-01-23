#!/usr/bin/env python3
import os

inp = """\
%mem=10GB
%nprocshared=16
%chk=Pre.chk
#p rhf/cc-pVTZ Pop=Hirshfeld

name

0 1
"""

charges = {
    'choline': '1',
    'c4mim': '1',
    'acetate': '-1',
    'ac': '-1',
    'mes': '-1',
    'dhp': '-1',
    'h2po4': '-1',
    'water': '0'
}


def get_files():
    for file in os.listdir('.'):
        if file.endswith('.xyz'):
            return file

def add_xyz(file, string):
    with open(file, "r") as f:
        for line in f.readlines()[2:]:
            string += line
    string += '\n' # gaussian needs a blank line at end of file
    return string

def make_cm5_input(file):
    jobfile = """\
#!/bin/bash
#PBS -l walltime=00:45:00
#PBS -l ncpus=16
#PBS -l mem=10GB
#PBS -l jobfs=10GB
#PBS -l software=g16
#PBS -l wd

module load gaussian/g16b01
g16 < cm5.inp > cm5.out 2>&1"""

    name = file[:-4]
    inputfile = inp.replace('name', name)
    inputfile = add_xyz(file, inputfile)
    with open('cm5.inp', 'w') as res:
        res.write(inputfile)
    with open('cm5.job', 'w') as job:
        job.write(jobfile)


parent = os.getcwd()
for path, _, files in os.walk('.'):
    for file in files:
        if file == 'equil.xyz':
            os.chdir(path)
            if not os.path.exists('cm5'):
                os.mkdir('cm5')
            os.chdir('cm5')
            os.system(f'cp ../{file} {file}')
            print(f'Making cm5 inputs in {os.getcwd()}')
            make_cm5_input(file)
            os.chdir(parent)
