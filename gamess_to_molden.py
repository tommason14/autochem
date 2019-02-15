#!/usr/bin/env python3

import sys
import re

if len(sys.argv) != 3:
    sys.exit('Syntax: gamess_to_molden logfile newfile')
logfile, newfile = sys.argv[1:]

def freq_data(file):
    """Parses GAMESS hessian calculation log file for the frequency data"""
    regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
    results = {'Modes': [], 'Frequencies': [], 'Reduced Mass': [], 'Intensities': []} # keys used as headers for csv
    with open(file, "r") as f:
        for line in f.readlines():
            if re.search(regex, line):
                mode, vib, symmetry, mass, intensity = line.split()
                mode = int(mode)
                vib, mass, intensity = map(float, (vib, mass, intensity))
                results['Modes'].append(mode)
                results['Frequencies'].append(vib)
                results['Reduced Mass'].append(mass)
                results['Intensities'].append(intensity)

    for key, value in results.items():
        results[key] = value[6:] #3N-6, with the 6 at the start = trans or rot modes.
    return results

def write_file(data, filename):
    string = "[Molden Format]\n"
    string += '[FREQ]\n'
    for i in data['Frequencies']:
        string += f'{i}\n'
    string += '[INT]\n'
    for i in data['Intensities']:
        string += f'{i}\n'
    with open(filename, "w") as f:
        f.write(string)
    

def main(log, new_file):
    data = freq_data(log)
    write_file(data, new_file)

 
main(logfile, newfile)