#!/usr/bin/env python3

import sys
import re
import numpy as np

if len(sys.argv) != 2:
    sys.exit('Syntax: gamess_to_molden.py logfile')
logfile = sys.argv[1]
newfile = logfile.rsplit('.')[0] + '.molden'

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

def calc_type(file):
    """ Returns the type of GAMESS calculation submitted """
    runtyp = None
    for line in read_file(file):
        if 'RUNTYP=OPTIMIZE' in line.upper():
            runtyp = "opt"
            break
        elif 'RUNTYP=ENERGY' in line.upper():
            runtyp = "spec"
            break
        if any(word in line.upper() for word in ('RUNTYP=HESSIAN', 'RUNTYPE=FMOHESS')):
            runtyp = "hessian"
            break
    return runtyp 
        

def find_init_coords(file):
    """
    Parses log file for initial coordinates.
    Output is different depending on whether FMO theory is used,
    and different methods for finding initial coordinates are required.
    """
    def get_atnum(symbol):

        ptable = {}
        #[symbol, mass, radius, connectors, vdw radii]
            #atomic weights from: http://www.ciaaw.org/atomic-weights.htm
        ptable[  0] = ['Xx',   0.00000, 0.00 ,  0, 0.430]
        ptable[  1] = [ 'H',   1.00798, 0.30 ,  1, 0.741]
        ptable[  2] = ['He',   4.00260, 0.99 ,  0, 0.880]
        ptable[  3] = ['Li',   6.96750, 1.52 ,  8, 0.550]
        ptable[  4] = ['Be',   9.01218, 1.12 ,  8, 1.030]
        ptable[  5] = [ 'B',  10.81350, 0.88 ,  6, 0.900]
        ptable[  6] = [ 'C',  12.01060, 0.77 ,  4, 0.880]
        ptable[  7] = [ 'N',  14.00685, 0.70 ,  3, 0.880]
        ptable[  8] = [ 'O',  15.99940, 0.66 ,  2, 0.840]
        ptable[  9] = [ 'F',  18.99840, 0.64 ,  1, 0.815]
        ptable[ 10] = ['Ne',  20.17970, 1.60 ,  0, 1.170]
        ptable[ 11] = ['Na',  22.98977, 1.86 ,  8, 1.300]
        ptable[ 12] = ['Mg',  24.30550, 1.60 ,  8, 1.550]
        ptable[ 13] = ['Al',  26.98154, 1.43 ,  8, 1.400]
        ptable[ 14] = ['Si',  28.08500, 1.17 ,  8, 1.250]
        ptable[ 15] = [ 'P',  30.97376, 1.10 ,  8, 1.220]
        ptable[ 16] = [ 'S',  32.06750, 1.04 ,  2, 1.190]
        ptable[ 17] = ['Cl',  35.45150, 0.99 ,  1, 0.995]
        ptable[ 18] = ['Ar',  39.94800, 1.92 ,  0, 1.530]
        ptable[ 19] = [ 'K',  39.09830, 2.31 ,  8, 1.190]
        ptable[ 20] = ['Ca',  40.07800, 1.97 ,  8, 1.640]
        ptable[ 21] = ['Sc',  44.95591, 1.60 ,  8, 1.670]
        ptable[ 22] = ['Ti',  47.86700, 1.46 ,  8, 1.530]
        ptable[ 23] = [ 'V',  50.94150, 1.31 ,  8, 1.550]
        ptable[ 24] = ['Cr',  51.99610, 1.25 ,  8, 1.555]
        ptable[ 25] = ['Mn',  54.93804, 1.29 ,  8, 1.540]
        ptable[ 26] = ['Fe',  55.84500, 1.26 ,  8, 1.530]
        ptable[ 27] = ['Co',  58.93319, 1.25 ,  8, 1.700]
        ptable[ 28] = ['Ni',  58.69340, 1.24 ,  8, 1.720]
        ptable[ 29] = ['Cu',  63.54600, 1.28 ,  8, 1.650]
        ptable[ 30] = ['Zn',  65.38000, 1.33 ,  8, 1.420]
        ptable[ 31] = ['Ga',  69.72300, 1.41 ,  8, 1.370]
        ptable[ 32] = ['Ge',  72.63000, 1.22 ,  8, 1.410]
        ptable[ 33] = ['As',  74.92159, 1.21 ,  8, 1.420]
        ptable[ 34] = ['Se',  78.97100, 1.17 ,  8, 1.410]
        ptable[ 35] = ['Br',  79.90400, 1.14 ,  1, 1.069]
        ptable[ 36] = ['Kr',  83.79800, 1.97 ,  0, 1.670]
        ptable[ 37] = ['Rb',  85.46780, 2.44 ,  8, 1.320]
        ptable[ 38] = ['Sr',  87.62000, 2.15 ,  8, 1.980]
        ptable[ 39] = [ 'Y',  88.90584, 1.80 ,  8, 1.760]
        ptable[ 40] = ['Zr',  91.22400, 1.57 ,  8, 1.680]
        ptable[ 41] = ['Nb',  92.90637, 1.41 ,  8, 1.670]
        ptable[ 42] = ['Mo',  95.95000, 1.36 ,  8, 1.550]
        ptable[ 43] = ['Tc',  98.00000, 1.35 ,  8, 1.600]
        ptable[ 44] = ['Ru', 101.07000, 1.33 ,  8, 1.650]
        ptable[ 45] = ['Rh', 102.90550, 1.34 ,  8, 1.700]
        ptable[ 46] = ['Pd', 106.42000, 1.38 ,  8, 1.790]
        ptable[ 47] = ['Ag', 107.86820, 1.44 ,  8, 1.890]
        ptable[ 48] = ['Cd', 112.41400, 1.49 ,  8, 1.830]
        ptable[ 49] = ['In', 114.81800, 1.66 ,  8, 1.660]
        ptable[ 50] = ['Sn', 118.71000, 1.62 ,  8, 0.000]
        ptable[ 51] = ['Sb', 121.76000, 1.41 ,  8, 0.000]
        ptable[ 52] = ['Te', 127.60000, 1.37 ,  8, 0.000]
        ptable[ 53] = [ 'I', 126.90447, 1.33 ,  1, 0.000]
        ptable[ 54] = ['Xe', 131.29300, 2.17 ,  0, 0.000]
        ptable[ 55] = ['Cs', 132.90545, 2.62 ,  8, 0.000]
        ptable[ 56] = ['Ba', 137.32700, 2.17 ,  8, 0.000]
        ptable[ 57] = ['La', 138.90547, 1.88 ,  8, 0.000]
        ptable[ 58] = ['Ce', 140.11600, 1.818,  8, 0.000]
        ptable[ 59] = ['Pr', 140.90766, 1.824,  8, 0.000]
        ptable[ 60] = ['Nd', 144.24200, 1.814,  8, 0.000]
        ptable[ 61] = ['Pm', 145.00000, 1.834,  8, 0.000]
        ptable[ 62] = ['Sm', 150.36000, 1.804,  8, 0.000]
        ptable[ 63] = ['Eu', 151.96400, 2.084,  8, 0.000]
        ptable[ 64] = ['Gd', 157.25000, 1.804,  8, 0.000]
        ptable[ 65] = ['Tb', 158.92535, 1.773,  8, 0.000]
        ptable[ 66] = ['Dy', 162.50000, 1.781,  8, 0.000]
        ptable[ 67] = ['Ho', 164.93033, 1.762,  8, 0.000]
        ptable[ 68] = ['Er', 167.25900, 1.761,  8, 0.000]
        ptable[ 69] = ['Tm', 168.93422, 1.759,  8, 0.000]
        ptable[ 70] = ['Yb', 173.04500, 1.922,  8, 0.000]
        ptable[ 71] = ['Lu', 174.96680, 1.738,  8, 0.000]
        ptable[ 72] = ['Hf', 178.49000, 1.57 ,  8, 0.000]
        ptable[ 73] = ['Ta', 180.94788, 1.43 ,  8, 0.000]
        ptable[ 74] = [ 'W', 183.84000, 1.37 ,  8, 0.000]
        ptable[ 75] = ['Re', 186.20700, 1.37 ,  8, 0.000]
        ptable[ 76] = ['Os', 190.23000, 1.34 ,  8, 0.000]
        ptable[ 77] = ['Ir', 192.21700, 1.35 ,  8, 0.000]
        ptable[ 78] = ['Pt', 195.08400, 1.38 ,  8, 0.000]
        ptable[ 79] = ['Au', 196.96657, 1.44 ,  8, 0.000]
        ptable[ 80] = ['Hg', 200.59200, 1.52 ,  8, 0.000]
        ptable[ 81] = ['Tl', 204.38350, 1.71 ,  8, 0.000]
        ptable[ 82] = ['Pb', 207.20000, 1.75 ,  8, 0.000]
        ptable[ 83] = ['Bi', 208.98040, 1.70 ,  8, 0.000]
        ptable[ 84] = ['Po', 209.00000, 1.40 ,  8, 0.000]
        ptable[ 85] = ['At', 210.00000, 1.40 ,  1, 0.000]
        ptable[ 86] = ['Rn', 222.00000, 2.40 ,  0, 0.000]
        ptable[ 87] = ['Fr', 223.00000, 2.70 ,  8, 0.000]
        ptable[ 88] = ['Ra', 226.00000, 2.20 ,  8, 0.000]
        ptable[ 89] = ['Ac', 227.00000, 2.00 ,  8, 0.000]
        ptable[ 90] = ['Th', 232.03770, 1.79 ,  8, 0.000]
        ptable[ 91] = ['Pa', 231.03588, 1.63 ,  8, 0.000]
        ptable[ 92] = [ 'U', 238.02891, 1.56 ,  8, 0.000]
        ptable[ 93] = ['Np', 237.00000, 1.55 ,  8, 0.000]
        ptable[ 94] = ['Pu', 244.00000, 1.59 ,  8, 0.000]
        ptable[ 95] = ['Am', 243.00000, 1.73 ,  8, 0.000]
        ptable[ 96] = ['Cm', 247.00000, 1.74 ,  8, 0.000]
        ptable[ 97] = ['Bk', 247.00000, 1.70 ,  8, 0.000]
        ptable[ 98] = ['Cf', 251.00000, 1.86 ,  8, 0.000]
        ptable[ 99] = ['Es', 252.00000, 1.86 ,  8, 0.000]
        ptable[100] = ['Fm', 257.00000, 2.00 ,  8, 0.000]
        ptable[101] = ['Md', 258.00000, 2.00 ,  8, 0.000]
        ptable[102] = ['No', 259.00000, 2.00 ,  8, 0.000]
        ptable[103] = ['Lr', 266.00000, 2.00 ,  8, 0.000]
        ptable[104] = ['Rf', 267.00000, 2.00 ,  8, 0.000]
        ptable[105] = ['Db', 268.00000, 2.00 ,  8, 0.000]
        ptable[106] = ['Sg', 269.00000, 2.00 ,  8, 0.000]
        ptable[107] = ['Bh', 270.00000, 2.00 ,  8, 0.000]
        ptable[108] = ['Hs', 277.00000, 2.00 ,  8, 0.000]
        ptable[109] = ['Mt', 278.00000, 2.00 ,  8, 0.000]
        ptable[110] = ['Ds', 281.00000, 2.00 ,  8, 0.000]
        ptable[111] = ['Rg', 282.00000, 2.00 ,  8, 0.000]
        ptable[112] = ['Cn', 285.00000, 2.00 ,  8, 0.000]
        ptable[113] = ['Nh', 286.00000, 2.00 ,  8, 0.000]
        ptable[114] = ['Fl', 289.00000, 2.00 ,  8, 0.000]
        ptable[115] = ['Mc', 290.00000, 2.00 ,  8, 0.000]
        ptable[116] = ['Lv', 293.00000, 2.00 ,  8, 0.000]
        ptable[117] = ['Ts', 294.00000, 2.00 ,  8, 0.000]
        ptable[118] = ['Og', 294.00000, 2.00 ,  8, 0.000]
    
        for k, v in ptable.items():
            if symbol == v[0]:
                return k

    bohrs = []
    angs = []

    regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"

    fmo = False
    for line in read_file(file):
        if 'FMO' in line:
            fmo = True
            break
    
    bohr_to_angs = 0.52917724899
    angs_to_bohr = 1.8897259886

    if fmo:
        calc = calc_type(file)
        if calc == "hessian":
            found = False
            # looking for this: 
            # 10  H                5.861571       6.324233       6.608315
            reg = '^\s*[0-9]{1,4}\s*[A-z]{1,2}(\s*-?[0-9]*.[0-9]*){3}$'
            atoms = []
            for line in read_file(file):
                if re.search('\s*COORD\s*0\s*NUCLEAR\s*COORDINATES', line):
                    found = True
                if found and 'TOTAL CPU TIME' in line:
                    break
                if found:
                    if re.search(reg, line):
                        _, sym, x, y, z = line.split()
                        x, y, z = map(float, (x, y, z))
                        x, y, z = map(lambda num: num * bohr_to_angs, (x, y, z))
                        atnum = get_atnum(sym)
                        atoms.append([sym, atnum, x, y, z])

            for index, atom in enumerate(atoms, 1):
                sym, atnum, x, y, z = atom
                angs.append(f"{sym:^3}{index:>4}{atnum:>4}{x:>12.6f}{y:>12.6f}{z:>12.6f}")
                x, y, z = map(lambda num : num * angs_to_bohr, (x, y, z))
                bohrs.append(f"{sym:^3}{x:>12.6f}{y:>12.6f}{z:>12.6f}")

        else:
            found = False
            # looking for this: 
            #C           6.0    -2.9137945280       -0.5434462158        0.0751487230
            reg = '^\s*[A-z]{1,2}(\s*-?[0-9]*.[0-9]*){4}$'
            atoms = []
            for line in read_file(file):
                if 'ATOM   CHARGE       X              Y              Z' in line:
                    found = True
                if found and 'TOTAL CPU TIME' in line:
                    break
                if found:
                    if re.search(reg, line):
                        sym, atnum, x, y, z = line.split()
                        x, y, z = map(float, (x, y, z))
                        atnum = get_atnum(sym)
                        atoms.append([sym, atnum, x, y, z])

            for index, atom in enumerate(atoms, 1):
                sym, atnum, x, y, z = atom
                angs.append(f"{sym:^3}{index:>4}{atnum:>4}{x:>12.6f}{y:>12.6f}{z:>12.6f}")
                x, y, z = map(lambda num : num * angs_to_bohr, (x, y, z))
                bohrs.append(f"{sym:^3}{x:>12.6f}{y:>12.6f}{z:>12.6f}")

    else:
        found = False
        # looking for this: 
        #C           6.0    -2.9137945280       -0.5434462158        0.0751487230
        reg = '^\s*[A-z]{1,2}(\s*-?[0-9]*.[0-9]*){4}$'
        atoms = []
        for line in read_file(file):
            if 'CHARGE         X                   Y                   Z' in line:
                found = True
            if found and line is '\n':
                break
            if found:
                if re.search(reg, line):
                    sym, atnum, x, y, z = line.split()
                    x, y, z = map(float, (x, y, z))
                    x, y, z = map(lambda num: num * bohr_to_angs, (x, y, z))
                    atnum = float(atnum)
                    atnum = int(atnum)
                    atoms.append([sym, atnum, x, y, z])

        for index, atom in enumerate(atoms, 1):
            sym, atnum, x, y, z = atom
            angs.append(f"{sym:^3}{index:>4}{atnum:>4}{x:>12.6f}{y:>12.6f}{z:>12.6f}")
            x, y, z = map(lambda num : num * angs_to_bohr, (x, y, z))
            bohrs.append(f"{sym:^3}{x:>12.6f}{y:>12.6f}{z:>12.6f}")

    return bohrs, angs

def find_geometries(file):
    """
    Parses GAMESS optimisations for geometries after every iteration. Returns a list of lists.
    """
    geoms = []
    atom_regex = '^\s*[A-z]+(\s*-?[0-9]+.[0-9]+){4}'
    found = False
    iteration = []
    scfs = []
    max_forces = []
    rms_forces = []
    for line in read_file(file):
        if 'EQUILIBRIUM GEOMETRY' in line:
            break # prints out coords of last iteration, not needed twice
        if 'YOU SHOULD RESTART' in line:
            break # prints out coords of last iteration, not needed twice
        if 'COORDINATES OF ALL ATOMS ARE' in line:
            found = True
            # reset atoms
        if found:
            if re.search(atom_regex, line):
                sym, atnum, x, y, z = line.split()
                x, y, z = map(float, (x, y, z))
                iteration.append(f"{sym:^3}{x:>12.6f}{y:>12.6f}{z:>12.6f}")
        if 'STEP CPU TIME' in line: #collected all atoms for that iteration
            found = False
        if not found: # finished that iteration
            if len(iteration) > 0:
                iteration = [f"  {len(iteration)}", ""] + iteration # add length and blank line
                for item in iteration:
                    geoms.append(item)
                iteration = []
        if 'NSERCH:' in line:
            line = line.split() 
            for ind, val in enumerate(line):
                if 'E=' in val:
                    scfs.append(line[ind + 1])
                elif 'MAX=' in val:
                    max_forces.append(line[ind + 1])
                elif 'R.M.S.=' in val:
                    rms_forces.append(line[ind + 1])
    
    options = {'energy': scfs,
               'max-force': max_forces,
               'rms-force': rms_forces}
    
    conv = []
    for k, v in options.items():
        conv.append(k)
        for val in v:
            conv.append(f"   {val}")

    return geoms, conv


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

def collect_into_dict(*,init_coords_bohr = None, init_coords_angs = None, wavenumbers = None,
intensities = None, vibrations = None, geometries = None, geom_convergence = None): 
    """
    Returns a dictionary whereby the keys are used by molden as a delimeter to define a new section.
    Values are lists of lines to write to the file, with no newline characters.
    Only necessary parameters are returned
    """
    options = {}
    options['Atoms'] = init_coords_angs
    options['GEOMETRIES'] = geometries
    options['GEOCONV'] = geom_convergence
    options['FREQ'] = wavenumbers
    options['INT'] = intensities
    options['FR-COORD'] = init_coords_bohr
    options['FR-NORM-COORD'] = vibrations

    ret = {}
    for k, v in options.items():
        if v is not None:
            ret[k] = v
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

def optimisation_params(file):
    """Finds parameters relevant to a GAMESS optimisation"""
    bohrs, angs = find_init_coords(file)
    num_atoms = len(angs)
    geometries, geom_conv = find_geometries(file)
    data = collect_into_dict(init_coords_bohr = bohrs, 
                             init_coords_angs = angs, 
                             geometries = geometries,
                             geom_convergence = geom_conv)
    return data


def hessian_params(file):
    """Finds parameters relevant to a GAMESS hessian calculation"""
    bohrs, angs = find_init_coords(log)
    num_atoms = len(angs)
    vibs, waves, ints = get_vibrations(num_atoms, log)
    data = collect_into_dict(init_coords_bohr = bohrs, 
                             init_coords_angs = angs, 
                             wavenumbers = waves, 
                             intensities = ints,
                             vibrations = vibs)
    return data

def main(log, new_file):

    calc = calc_type(log)
    print(calc)
    if calc == "opt":
        data = optimisation_params(log)
        write_file(data, new_file)
    elif calc == "hessian":
        data = hessian_params(log)
        write_file(data, new_file)

 
main(logfile, newfile)
