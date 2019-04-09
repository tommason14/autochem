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
            for line in f:
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
    get_atnum = {} # dict lookup faster than func. call
    get_atnum['Xx'] =   0
    get_atnum[ 'H'] =   1
    get_atnum['He'] =   2
    get_atnum['Li'] =   3
    get_atnum['Be'] =   4
    get_atnum[ 'B'] =   5
    get_atnum[ 'C'] =   6
    get_atnum[ 'N'] =   7
    get_atnum[ 'O'] =   8
    get_atnum[ 'F'] =   9
    get_atnum['Ne'] =  10
    get_atnum['Na'] =  11
    get_atnum['Mg'] =  12
    get_atnum['Al'] =  13
    get_atnum['Si'] =  14
    get_atnum[ 'P'] =  15
    get_atnum[ 'S'] =  16
    get_atnum['Cl'] =  17
    get_atnum['Ar'] =  18
    get_atnum[ 'K'] =  19
    get_atnum['Ca'] =  20
    get_atnum['Sc'] =  21
    get_atnum['Ti'] =  22
    get_atnum[ 'V'] =  23
    get_atnum['Cr'] =  24
    get_atnum['Mn'] =  25
    get_atnum['Fe'] =  26
    get_atnum['Co'] =  27
    get_atnum['Ni'] =  28
    get_atnum['Cu'] =  29
    get_atnum['Zn'] =  30
    get_atnum['Ga'] =  31
    get_atnum['Ge'] =  32
    get_atnum['As'] =  33
    get_atnum['Se'] =  34
    get_atnum['Br'] =  35
    get_atnum['Kr'] =  36
    get_atnum['Rb'] =  37
    get_atnum['Sr'] =  38
    get_atnum[ 'Y'] =  39
    get_atnum['Zr'] =  40
    get_atnum['Nb'] =  41
    get_atnum['Mo'] =  42
    get_atnum['Tc'] =  43
    get_atnum['Ru'] =  44
    get_atnum['Rh'] =  45
    get_atnum['Pd'] =  46
    get_atnum['Ag'] =  47
    get_atnum['Cd'] =  48
    get_atnum['In'] =  49
    get_atnum['Sn'] =  50
    get_atnum['Sb'] =  51
    get_atnum['Te'] =  52
    get_atnum[ 'I'] =  53
    get_atnum['Xe'] =  54
    get_atnum['Cs'] =  55
    get_atnum['Ba'] =  56
    get_atnum['La'] =  57
    get_atnum['Ce'] =  58
    get_atnum['Pr'] =  59
    get_atnum['Nd'] =  60
    get_atnum['Pm'] =  61
    get_atnum['Sm'] =  62
    get_atnum['Eu'] =  63
    get_atnum['Gd'] =  64
    get_atnum['Tb'] =  65
    get_atnum['Dy'] =  66
    get_atnum['Ho'] =  67
    get_atnum['Er'] =  68
    get_atnum['Tm'] =  69
    get_atnum['Yb'] =  70
    get_atnum['Lu'] =  71
    get_atnum['Hf'] =  72
    get_atnum['Ta'] =  73
    get_atnum[ 'W'] =  74
    get_atnum['Re'] =  75
    get_atnum['Os'] =  76
    get_atnum['Ir'] =  77
    get_atnum['Pt'] =  78
    get_atnum['Au'] =  79
    get_atnum['Hg'] =  80
    get_atnum['Tl'] =  81
    get_atnum['Pb'] =  82
    get_atnum['Bi'] =  83
    get_atnum['Po'] =  84
    get_atnum['At'] =  85
    get_atnum['Rn'] =  86
    get_atnum['Fr'] =  87
    get_atnum['Ra'] =  88
    get_atnum['Ac'] =  89
    get_atnum['Th'] =  90
    get_atnum['Pa'] =  91
    get_atnum[ 'U'] =  92
    get_atnum['Np'] =  93
    get_atnum['Pu'] =  94
    get_atnum['Am'] =  95
    get_atnum['Cm'] =  96
    get_atnum['Bk'] =  97
    get_atnum['Cf'] =  98
    get_atnum['Es'] =  99
    get_atnum['Fm'] = 100
    get_atnum['Md'] = 101
    get_atnum['No'] = 102
    get_atnum['Lr'] = 103
    get_atnum['Rf'] = 104
    get_atnum['Db'] = 105
    get_atnum['Sg'] = 106
    get_atnum['Bh'] = 107
    get_atnum['Hs'] = 108
    get_atnum['Mt'] = 109
    get_atnum['Ds'] = 110
    get_atnum['Rg'] = 111
    get_atnum['Cn'] = 112
    get_atnum['Nh'] = 113
    get_atnum['Fl'] = 114
    get_atnum['Mc'] = 115
    get_atnum['Lv'] = 116
    get_atnum['Ts'] = 117
    get_atnum['Og'] = 118
    
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
                        atnum = get_atnum[sym]
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
                        atnum = get_atnum[sym]
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

def find_geometries(file, num_atoms):
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
       
    
    # trim output so the scf convergence doesn't show a spike
    # to 0 (happends if len(scfs) is not equal to geoms
    # - depends on output from GAMESS, sometimes prints
    # extra coords

    num_iterations = int(len(geoms) / (num_atoms + 2)) 
    # + 2 because add num_atoms and blank line
    # for every iteration
    if num_iterations > len(scfs):
        geoms = geoms[:len(scfs) * (num_atoms + 2)] # trim to fit

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
    geometries, geom_conv = find_geometries(file, num_atoms)
    data = collect_into_dict(init_coords_bohr = bohrs, 
                             init_coords_angs = angs, 
                             geometries = geometries,
                             geom_convergence = geom_conv)
    return data


def hessian_params(file):
    """Finds parameters relevant to a GAMESS hessian calculation"""
    bohrs, angs = find_init_coords(file)
    num_atoms = len(angs)
    vibs, waves, ints = get_vibrations(num_atoms, file)
    data = collect_into_dict(init_coords_bohr = bohrs, 
                             init_coords_angs = angs, 
                             wavenumbers = waves, 
                             intensities = ints,
                             vibrations = vibs)
    return data

def main(log, new_file):

    calc = calc_type(log)
    if calc == "opt":
        data = optimisation_params(log)
        write_file(data, new_file)
    elif calc == "hessian":
        data = hessian_params(log)
        write_file(data, new_file)

 
main(logfile, newfile)
