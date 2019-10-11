from .atom import Atom
from .utils import read_file, write_csv_from_dict
import os
import re
import subprocess
import sys

__all__ = ['thermo_data']

def get_filetype(file):
    """
    Returns either 'gamess' or 'gauss' depending on filetype of log file
    """
    for line in read_file(file):
        if 'gamess' in line.lower():
            return 'gamess'
        if 'gaussian' in line.lower():
            return 'gauss'   

def write_geom_input(atoms):
    """
    Writes 'geom.input' from the list of |Atom| objects passed in.
    """
    with open('geom.input', 'w') as new:
        for atom in atoms:
            new.write(f"{atom.symbol:5s} {int(atom.atnum):3} {atom.x:>15.10f} {atom.y:>15.10f} {atom.z:>15.10f} \n")

def thermo_initial_geom_gamess(file):
    """Parses GAMESS hessians for the initial geometry"""
    atoms = []
    regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"
    inp = file[:-3] + 'inp'
    for line in read_file(inp):
        if re.search(regex, line):
            sym, _, x, y, z = line.split()
            x, y, z = map(float, (x, y, z))
            atoms.append(Atom(symbol = sym, coords = (x, y, z)))
    write_geom_input(atoms)

def gauss_num_atoms(file):
    """
    Returns number of atoms in Gaussian calc
    """
    for line in read_file(file):
        if 'NAtoms=' in line:
            return int(line.split()[1])

def thermo_initial_geom_gauss(file):
    """
    Parses Gaussian frequency calculation log file for the initial geometry. Note that
    coordinates here are stored in .job files by default.
    Only works with xyz coordinates, not z-matrices.
    Also accounts for the fact that an optimisation might occur before the frequency job,
    so if 'Input orientation:' is found and atoms is not empty, the original atoms list is deleted."""
    atoms = []
    num_atoms = gauss_num_atoms(file)
    regex = "(\s+[0-9]+){3}(\s+-?[0-9]+\.[0-9]+){3}"
    found_freq = False
    found_coords = False
    for line in read_file(file):
        # deal with opt-freq or just freq
        if re.search('Freq$', line) or \
            'freq' in line.lower() and \
            'opt' not in line.lower():
            found_freq = True
        if 'Input orientation:' in line:
            found_coords = True
        if 'Distance matrix' in line:
            found_coords = False
        if found_freq and found_coords and re.search(regex, line):
            _, atnum, _, x, y, z = line.split()
            atnum, x, y, z = map(float, (atnum, x, y, z))
            atoms.append(Atom(atnum = atnum, coords = (x, y, z)))
    write_geom_input(atoms)

def freq_data_gamess(file, write_freqs_to_file = False):
    """Parses GAMESS hessian log files for frequency data"""
    regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
    found_region = False
    modes = []
    freqs = []
    ints  = []
    with open(file, "r") as f:
        for line in f:
            if 'MODE FREQ(CM**-1)  SYMMETRY  RED. MASS  IR INTENS.' in line:
                found_region = True
            if line is '\n':
                found_region = False
            if found_region:
                if re.search(regex, line):
                    mode, vib, *_,  intensity = line.split()
                    mode = int(mode)
                    vib, intensity = map(float, (vib, intensity))
                    modes.append(mode)
                    freqs.append(vib)
                    ints.append(intensity)

    results = {'Modes'                          : modes, 
               'Frequencies [cm-1]'             : freqs, 
               'Intensities [Debye^2/(amu Å^2)]': ints} # keys used as headers for csv

    for key, value in results.items():
        results[key] = value[6:] #3N-6, with the 6 at the start = trans or rot modes.

    if write_freqs_to_file:
        with open("freq.out", "w") as output:
            for i in results['Frequencies [cm-1]']:
                output.write(f"{i:.3f}\n")
    return results
    
def freq_data_gauss(file, write_freqs_to_file = False):
    """Parses Gaussian frequency log files for frequency data"""
    # regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
    found_region = False
    freqs = []
    ints  = []
    with open(file, "r") as f:
        for line in f:
            if 'Frequencies --' in line:
                freqs += line.split()[2:]
            if 'IR Inten    --' in line:
                ints += line.split()[3:]
 
    freqs = [float(i) for i in freqs]
    ints = [float(i) for i in ints]
    modes = [i for i in range(len(freqs))]
    results = {'Modes'                          : modes, 
               'Frequencies [cm-1]'             : freqs, 
               'Intensities [Debye^2/(amu Å^2)]': ints} # keys used as headers for csv

    for key, value in results.items():
        results[key] = value[6:] #3N-6, with the 6 at the start = trans or rot modes.
    
    if write_freqs_to_file:
        with open("freq.out", "w") as output:
            for i in results['Frequencies [cm-1]']:
                output.write(f"{i:.3f}\n")
    return results



def run(file, mult, temp):
    """
    Calls thermo.exe with geom.input and freq.out written to the same directory.
    """
    thermo_exe = os.path.join(os.path.dirname(os.path.realpath(__file__)), f'thermo.exe')
    p = subprocess.Popen(thermo_exe, shell=True, 
    stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True) 
    newline = os.linesep
    commands = ['y', 'y', 'y', mult, temp]
    p.communicate(newline.join(commands))

def read_fort():
    with open('fort.10', 'r') as f:
        fort = [line for line in f.readlines()]
    return fort

def grep_data(fort):
    """Collects desired data from fort.10"""
    data = {}
    lookup = \
    {
        'TC h'   : 'TC',
        'S elec' : 'S elec',
        'S trans' : 'S trans',
        'S rot'  : 'S rot',
        'S vib'  : 'S vib',
        'Stot'   : 'S tot',
        'TC-'    : 'TC - TS'
    }
    for line in fort:
        if re.search('ZPVE.*=.*kJ', line):
            data['ZPVE'] = line.split()[2]
        for k, v in lookup.items():
            if k in line:
                data[v] = line.split()[-1]
    return data

def cleanup():
    os.system('rm fort.10 moments geom.input freq.out')

def setup_and_run_fortran_script(file, mult, temp):
    """
    Runs fortran script to produce 'fort.10' files etc...
    """
    filetype=get_filetype(file)
    if filetype == 'gamess': 
        thermo_initial_geom_gamess(file)
        freq_data_gamess(file, write_freqs_to_file = True)
        run(file, mult, temp)
    if filetype == 'gauss':
        thermo_initial_geom_gauss(file)
        freq_data_gauss(file, write_freqs_to_file = True)
        run(file, mult, temp)

def thermo_data(file, mult, temp):
    """
    Uses a fortran script to produce thermochemical data for GAMESS Hessian calculations and
    GAUSSIAN frequency calculations- the results produced in the log file have been shown to be 
    inaccurate. At < 300 cm⁻¹, rigid rotor fails, and the fortran code implements hindered rotor.
    """
    setup_and_run_fortran_script(file, mult, temp)
    fort = read_fort()    
    data = grep_data(fort)
    cleanup()
    return data

