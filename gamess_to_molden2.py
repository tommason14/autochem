#!/usr/bin/env python3


import sys
import re
import numpy as np

if len(sys.argv) != 3:
    sys.exit('Syntax: gamess_to_molden logfile newfile')
logfile, newfile = sys.argv[1:]


def read_file(file):
    """
    Returns an iterator over a file object
    """
    
    with open(file, "r") as f:
        try:
            for line in f.readlines():
                yield line
        except UnicodeDecodeError:
            pass

def find_init_coords(file):
    """
    Parses input file for initial coordinates
    """

    bohrs = []
    angs = []

    inp_file = file[:-3] + 'inp'

    regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"

    init_coords = []
    for line in read_file(inp_file):
       if re.search(regex, line): 
            sym, atnum, x, y, z = line.split()
            atnum = float(atnum)
            x, y, z = map(float, (x, y, z))
            init_coords.append([sym, atnum, x, y, z])

    def angstrom_to_bohr(num):
        return num * 1.8897259886
    
    for index, atom in enumerate(init_coords, 1):
        sym, atnum, x, y, z = atom
        atnum = int(atnum)
        angs.append(f"{sym:^3}{index:>4}{atnum:>4}{x:>12.6f}{y:>12.6f}{z:>12.6f}")
        x, y, z = map(angstrom_to_bohr, (x, y, z))
        bohrs.append(f"{sym:^3}{x:>12.6f}{y:>12.6f}{z:>12.6f}")
                
    return bohrs, angs

def freq_data(file):
    """
    Parses GAMESS hessian calculation log file for the frequency data
    """

    regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
    results = {'Modes': [], 'Frequencies': [], 'Reduced Mass': [], 'Intensities': []}
    for line in read_file(file):
        if re.search(regex, line):
            mode, vib, symmetry, mass, intensity = line.split()
            mode = int(mode)
            vib, mass, intensity = map(float, (vib, mass, intensity))
            results['Modes'].append(mode)
            results['Frequencies'].append(vib)
            results['Reduced Mass'].append(mass)
            results['Intensities'].append(intensity)
    return results

def find_normal_coords_approx(file):
    """
    Parses log for the changes in atom positions that occurs with each vibration.
    Regex for lines matches in multiple places, so just locating the vibrations within a
    section of the file- see `refine selection` for how to get the actual vibrations
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

def refine_selection(lst):
    """
    Remove unnecessary data from a list of lines from log file
    """

    # frequencies = []
    # intensities  = []
    refined = []
    xvib_regex = '^\s*?[0-9]*\s*[A-z]{1,2}\s*[X]' 

    # remove unnecessary lines
    found_vibs = False
    for line in lst:
        # if re.search('\s*FREQUENCY:', line):
            # frequencies += line.split()[1:]
        # if re.search('\s*IR\sINTENSITY:', line):
            # intensities += line.split()[2:]
        if re.search(xvib_regex, line):
            found_vibs = True
        if '\n' is line:
            found_vibs = False
        if found_vibs:
            refined.append(line.strip()) # rm \n
    # return refined, frequencies, intensities
    # works, but refactoring needed.
    # for now, settle for opening file one more time

    return refined
    

def find_normal_coords(num_atoms, file):
    """
    Parses log for the changes in atom positions that occurs with each vibration
    """

    modes, approx = find_normal_coords_approx(file)
    vibs = {}
    xvib_regex = '^\s*?[0-9]*\s*[A-z]{1,2}\s*[X]' 
    xvibs_per_line = []
    yvibs_per_line = []
    zvibs_per_line = []

    refined = refine_selection(approx)

    for line in refined:
        if re.search(xvib_regex, line):
            xvibs_per_line.append(line.split()[3:])
        if line.startswith('Y'):
            yvibs_per_line.append(line.split()[1:])
        if line.startswith('Z'):
            zvibs_per_line.append(line.split()[1:])

    for a in range(1, num_atoms + 1):
        vibs[a] = {'x': [], 'y' : [], 'z': []}
 
    #looping over list of lists- n x number of atoms
    #
    #     [ [ ] , [ ] , [ ] , [ ] , [ ] , [ ] ]  
    #        ^     ^     ^     ^     ^     ^
    # atom:  1     2     3     1     2     3
    
    # changes in positions of atoms
    def add_deltas(d, lst, dim):
        """Add vibrations to the correct atom and correct dimension"""
        atom_count = 0
        for num, line in enumerate(lst):
            atom_count += 1
            d[atom_count][dim] += line
            # reset count when reached max number of atoms
            if atom_count == num_atoms:
                atom_count = 0 
        return d


    vibs = add_deltas(vibs, xvibs_per_line, 'x')
    vibs = add_deltas(vibs, yvibs_per_line, 'y')
    vibs = add_deltas(vibs, zvibs_per_line, 'z')

    return vibs

   # {atom 1: {'x': [mode1, mode2, mode3...],
   #           'y': [mode1, mode2, mode3...],
   #           'z': [mode1, mode2, mode3...]
   #          },
   #  atom 2: {'x': [mode1, mode2, mode3...],
   #           'y': [mode1, mode2, mode3...],
   #           'z': [mode1, mode2, mode3...],
   #          }...


def tidy_normal_modes(data):
    """
    Rearrange the normal modes into a desirable format
    """
   
    num_atoms = len(data)
    num_modes = len(data[1]['x']) # look at any list, all the same length
    num_dimensions = 3

    tidy = np.zeros((num_modes, num_atoms, num_dimensions))

    for a, atom in enumerate(data.keys()):
        for d, dim in enumerate(data[atom].keys()):
            for mode, coord in enumerate(data[atom][dim]):
                tidy[mode,a, d] = coord
   
    return tidy

def pretty_print_vibs(vibrations):
    """
    Prints vibrations in a molden-compatible format
    """
    vibs = []
    for vib, v in enumerate(vibrations):
        vibs.append(f'vibration     {vib + 1}') 
        for atom in v:
            x, y, z = atom
            vibs.append(f"{x:>12.6f}{y:>13.6f}{z:>13.6f}")

def get_vibrations(num_atoms, log):
    """
    Control flow for finding vibrations.
    """
    normals = find_normal_coords(num_atoms, log)
    vibs = tidy_normal_modes(normals)
    vibs = pretty_print_vibs(vibs)
    return vibs

def collect_into_dict(*,init_coords_bohr = None, init_coords_angs = None, freq_data = None, vibrations = None): 
    """
    Returns a dictionary whereby the keys are used by molden as a delimeter to define a new section.
    Values are lists of lines    
    """

    ret = {}
    ret['Atoms'] = init_coords_angs
    ret['FREQ'] = freq_data['Frequencies']
    ret['INT'] = freq_data['Intensities']
    ret['FR-COORD'] = init_coords_bohr
    ret['FR-NORM-COORD'] = vibs
    return ret


def write_file(data, filename):
    """
    Writes a molden-readable file
    """

    lst = ["[Molden Format]\n"]

    for header, lines in data.items():
        lst.append(f'[{header}]\n')
        for line in lines:
            lst.append(f'{line}\n')

    with open(filename, "w") as f:
        for line in lst:
            f.write(line)

def main(log, new_file):

    bohrs, angs = find_init_coords(log)
    num_atoms = len(angs)
    freq = freq_data(log)
    vibs = get_vibrations(num_atoms, log)
    data = collect_into_dict(init_coords_bohr = bohrs, 
                             init_coords_angs = angs, 
                             freq_data = freq, 
                             vibrations = vibs)
    write_file(data, new_file)

 
main(logfile, newfile)
