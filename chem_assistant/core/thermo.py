from .atom import Atom

import os
import re
import subprocess

__all__ = ['thermo_data']

def thermo_initial_geom(file):
    atoms = []
    regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"
    found = False
    with open(file, "r") as f:
        for line in f.readlines():
            if 'CHARGE         X                   Y                   Z' in line:
                found = True
            if found:
                if re.search(regex, line):
                    sym, _, x, y, z = line.split()
                    x, y, z = map(float, (x, y, z))
                    atoms.append(Atom(symbol = sym, coords = (x, y, z)))
            if line is '\n':
                found = False
    with open('geom.input', 'w') as new:
        for atom in atoms:
            new.write(f"{atom.symbol:5s} {str(atom.atnum):3s} {atom.x:>15.10f} {atom.y:>15.10f} {atom.z:>15.10f} \n")

def thermo_freqs(file):
    regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z]{1,9}\s*[0-9]{1,3}\.[0-9]{1,9}\s*[0-9]{1,9}\.[0-9]{1,9}$'
    vibs = []
    with open(file, "r") as f:
        for line in f.readlines():
            if re.search(regex, line):
                vibs.append(line.split()[1])
    vibs = vibs[6:] #3N-6, with the 6 at the start = trans/rot.
    with open("freq.out", "w") as output:
        for i in vibs:
            output.write(i + "\n")

def run(file):
    p = subprocess.Popen(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'thermo-gamess.exe'), shell=True, stdin=subprocess.PIPE,
    stdout=subprocess.PIPE, universal_newlines=True)
    newline = os.linesep # [1]
    commands = ['y', 'y', 'y', '1', '298.15']
    p.communicate(newline.join(commands))
    print("\033[30;46mThermodynamics for", file +"\033[0m")
    os.system("cat fort.10")

def cleanup():
    os.system('rm fort.10 moments geom.input freq.out')

def thermo_data(file):
    """Uses a fortran script to produce thermochemical data for GAMESS Hessian calculations- the
results produced in the log file have been shown to be inaccurate."""
    thermo_initial_geom(file)
    thermo_freqs(file)
    run(file)
    cleanup()
