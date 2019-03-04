#!/usr/bin/env python3

import sys
import re
import numpy as np
import time


def timeit(f):

    def timed(*args, **kw):

        ts = time.time()
        result = f(*args, **kw)
        te = time.time()

        print('func:%r took: %2.4f sec' % \
          (f.__name__, te-ts))
        return result

    return timed


if len(sys.argv) != 2:
    sys.exit('Syntax: gamess_to_molden.py logfile')
logfile = sys.argv[1]
newfile = logfile[:-3] + 'molden'

def read_file(file):
    """
    Returns an iterator over a file object.
    """
    
    with open(file, "r") as f:
        try:
            for line in f.readlines():
                yield line
        except UnicodeDecodeError:
            pass

def find_init_coords(file):
    """
    Parses input file for initial coordinates.
    """

    bohrs = []
    angs = []

    inp_file = file[:-3] + 'inp'

    regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"

    init_coords = []
    try:
        for line in read_file(inp_file):
           if re.search(regex, line): 
                sym, atnum, x, y, z = line.split()
                atnum = float(atnum)
                x, y, z = map(float, (x, y, z))
                init_coords.append([sym, atnum, x, y, z])
    except FileNotFoundError:
        print(f'Error: Needs an input file in the same directory, {inp_file}')
        sys.exit()

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
    Parses GAMESS hessian calculation log file for the frequency data.
    Note: increased run time to more than 5 seconds! Now deprecated in favour of
    additional output from `refine_selection`.
    """
    regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
    results = {'Modes': [], 'Frequencies': [], 'Reduced Mass': [], 'Intensities': []}
    for line in read_file(file):
        if re.search(regex, line):
            mode, vib, symmetry, mass, intensity = line.split()
            mode = int(mode)
            vib, mass, intensity = map(float, (vib, mass, intensity))
            results['Frequencies'].append(vib)
            results['Intensities'].append(intensity)
    return results

def find_normal_coords_approx(file):
    """
    Parses log for the changes in atom positions that occurs with each vibration.
    Regex for lines matches in multiple places, so just locating the vibrations within a
    section of the file- see `refine selection` for how to get the actual vibrations.
    """

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
    Remove unnecessary data from a list of lines from log file.
    """

    frequencies = []
    intensities  = []
    refined = []
    xvib_regex = '^\s*?[0-9]*\s*[A-z]{1,2}\s*[X]' 

    found_vibs = False
    for line in lst:
        if re.search('\s*FREQUENCY:', line):
            frequencies += [freq for freq in line.split()[1:] if freq != 'I']
        if re.search('\s*IR\sINTENSITY:', line):
            intensities += line.split()[2:]
        if re.search(xvib_regex, line):
            found_vibs = True
        if '\n' is line:
            found_vibs = False
        if found_vibs:
            refined.append(line.strip()) # rm \n
    return refined, frequencies, intensities
    

def find_normal_coords(num_atoms, file):
    """
    Parses log for the changes in atom positions that occurs with each vibration.

    Example output, where the columns represent changes in position for each mode,
    and rows represent the X,Y and Z components for each atom in the system.

    1   C            X -0.00738929  0.15398404 -0.02729367 -0.08048432  0.09934430
                     Y  0.00055975  0.03053089  0.15133284  0.08248114  0.08159716
                     Z  0.00514209 -0.10069700  0.00686349  0.06067365  0.02132913
    """

    modes, approx = find_normal_coords_approx(file)
    vibs = {}
    xvib_regex = '^\s*?[0-9]*\s*[A-z]{1,2}\s*[X]' 
    xvibs_per_line = []
    yvibs_per_line = []
    zvibs_per_line = []

    refined, wavenumbers, intensities = refine_selection(approx)

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

    return vibs, wavenumbers, intensities

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
    Rearrange the normal modes into a desirable format.
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
    Prints vibrations in a molden-compatible format.
    """
    vibs = []
    for vib, v in enumerate(vibrations):
        vibs.append(f'vibration     {vib + 1}') 
        for atom in v:
            x, y, z = atom
            vibs.append(f"{x:>12.6f}{y:>13.6f}{z:>13.6f}")
    return vibs

def get_vibrations(num_atoms, log):
    """
    Control flow for finding vibrations.
    """
    normals, wavenumbers, intensities = find_normal_coords(num_atoms, log)
    tidy = tidy_normal_modes(normals)
    vibs = pretty_print_vibs(tidy)
    return vibs, wavenumbers, intensities

def collect_into_dict(*,init_coords_bohr = None, init_coords_angs = None, wavenumbers = None, intensities = None, vibrations = None): 
    """
    Returns a dictionary whereby the keys are used by molden as a delimeter to define a new section.
    Values are lists of lines to write to the file, with no newline characters.
    """

    ret = {}
    ret['Atoms'] = init_coords_angs
    ret['FREQ'] = wavenumbers
    ret['INT'] = intensities
    ret['FR-COORD'] = init_coords_bohr
    ret['FR-NORM-COORD'] = vibrations
    return ret

def write_file(data, filename):
    """
    Writes a molden-readable file.
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
    vibs, waves, ints = get_vibrations(num_atoms, log)
    data = collect_into_dict(init_coords_bohr = bohrs, 
                             init_coords_angs = angs, 
                             wavenumbers = waves, 
                             intensities = ints,
                             vibrations = vibs)
    write_file(data, new_file)

 
main(logfile, newfile)
