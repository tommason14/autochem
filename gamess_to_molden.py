#!/usr/bin/env python3

# imports {{{1

import sys
import re

if len(sys.argv) != 3:
    sys.exit('Syntax: gamess_to_molden logfile newfile')
logfile, newfile = sys.argv[1:]


# functions {{{1

def read_file(file):
    with open(file, "r") as f:
        try:
            for line in f.readlines():
                yield line
        except UnicodeDecodeError:
            pass

def find_init_coords(file):
    """Parses log for initial coordinates"""
    coords = []
    found = False
    regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"
    for line in read_file(file):
        if 'CHARGE         X                   Y                   Z' in line:
            found = True
        if found:
            if re.search(regex, line):
                sym, _, x, y, z = line.split()
                coords.append('  '.join([sym, x, y, z]))
        if line is '\n':
            found = False
    return coords

def freq_data(file):
    """Parses GAMESS hessian calculation log file for the frequency data"""
    regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
    results = {'Modes': [], 'Frequencies': [], 'Reduced Mass': [], 'Intensities': []} # keys used as headers for csv
    for line in read_file(file):
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

def find_normal_coords_approx(file):
    """
    Parses log for the changes in atom positions that occurs with each vibration.
    Regex for lines matches in multiple places, so just locating the vibrations within a
    section of the file
    """
    # Need to match multiple lines, so can't read line-by-line?
    normals = []
    num_modes = ''
    found = False
    
    for line in read_file(file):
        if 'ANALYZING SYMMETRY OF NORMAL MODES' in line:
            found = True
        if found:
            normals.append(line)
        if 'REFERENCE ON SAYVETZ CONDITIONS' in line:
            found = False
        if 'VIBRATIONAL MODES ARE USED IN THERMOCHEMISTRY' in line:
            num_modes = line.split()[2]
    return num_modes, normals


def find_normal_coords(file):
    """Parses log for the changes in atom positions that occurs with each vibration"""
    modes, approx = find_normal_coords_approx(file)
    refined = []
    for line in approx:
        if line.startswith(' '*29):
            pass # mode numbers         

    """
    Distinguish between lines starting with:
    -          1           2           3           4           5 (mode number)
    - FREQUENCY
    - SYMMETRY
    - REDUCED MASS
    - IR INTENSITY
    - Blank lines (in loop, pass)
    -  1   O            X 
    -                   Y (by regex \s{18}(Y|Z)
    -                   Z 
    -  2   H            X 
    -                   Y 
    -                   Z 
    -  3   H            X 
    -                   Y 
    -                   Z 
    """ 

    ### mode numbers ###
    




# collate {{{1 
def collect_into_dict(*,init_coords = None, freq_data = None): 
    """
    Returns a dictionary where the key is used by molden as a delimeter.
    Values are lists of lines    
    """
    ret = {}
    ret['FR-COORD'] = init_coords
    ret['FREQ'] = freq_data['Frequencies']
    ret['INT'] = freq_data['Intensities']
    return ret

# write to file {{{1

def write_file(data, filename):
    lst = ["[Molden Format]\n"]

    for header, lines in data.items():
        lst.append(f'[{header}]\n')
        for line in lines:
            lst.append(f'{line}\n')

    # string =     string += '[FREQ]\n'
    # for i in data['Frequencies']:
    #     string += f'{i}\n'
    # string += '[INT]\n'
    # for i in data['Intensities']:
    #     string += f'{i}\n'
    with open(filename, "w") as f:
        for line in lst:
            f.write(line)

# main {{{1
def main(log, new_file):
    init = find_init_coords(log)
    freq = freq_data(log)
    normals = find_normal_coords(log)
    data = collect_into_dict(init_coords = init, freq_data = freq)
    # for k, v in data.items():
    #     print(k)    
    #     for i in v: # elements of list
    #         print(i)
    # write_file(data, new_file)

 
main(logfile, newfile)
