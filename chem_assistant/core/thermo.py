from .atom import Atom
from .utils import write_csv_from_dict
import os
import re
import subprocess

__all__ = ['thermo_data', 'make_ir_spectra']

def thermo_initial_geom(file):
    """Parses GAMESS hessian calculation log file for the initial geometry"""
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

def freq_data(file, write_freqs_to_file = False):
    """Parses GAMESS hessian calculation log file for the frequency data"""
    regex = '[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$'
    results = {'Modes': [], 'Frequencies [cm-1]': [], 'Reduced Mass [amu]': [], 'Intensities [Debye^2/(amu Å^2)]': []} # keys used as headers for csv
    with open(file, "r") as f:
        for line in f.readlines():
            if re.search(regex, line):
                mode, vib, symmetry, mass, intensity = line.split()
                mode = int(mode)
                vib, mass, intensity = map(float, (vib, mass, intensity))
                results['Modes'].append(mode)
                results['Frequencies [cm-1]'].append(vib)
                results['Reduced Mass [amu]'].append(mass)
                results['Intensities [Debye^2/(amu Å^2)]'].append(intensity)

    for key, value in results.items():
        results[key] = value[6:] #3N-6, with the 6 at the start = trans or rot modes.

    if write_freqs_to_file:
        with open("freq.out", "w") as output:
            for i in results['vibs']:
                output.write(f"{i}\n")
    return results
    
def run(file):
    p = subprocess.Popen(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'thermo-gamess.exe'), shell=True, stdin=subprocess.PIPE,
    stdout=subprocess.PIPE, universal_newlines=True) 
    # os.path.dirname(os.path.realpath(__file__)) = dir of this file, irrespective of where the
    # python script is executed
    newline = os.linesep
    commands = ['y', 'y', 'y', '1', '298.15']
    p.communicate(newline.join(commands))
    # print("\033[30;46mThermodynamics for", file +"\033[0m")
    # os.system("cat fort.10")

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
        'S tran' : 'S tran',
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

def thermo_data(file):
    """Uses a fortran script to produce thermochemical data for GAMESS Hessian calculations- the
results produced in the log file have been shown to be inaccurate."""
    thermo_initial_geom(file)
    freq_data(file, write_freqs_to_file = True)
    run(file)
    fort = read_fort()    
    data = grep_data(fort)
    cleanup()
    return data

def make_ir_spectra(file):
    """Plot of wavenumber against intensities for vibrations found by diagonalisation of a computed
hessian matrix"""

    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt
    import seaborn as sns
    # Slow imports, so only import if needed

    res = freq_data(file)
    sns.set_style('darkgrid')
    sns.lineplot(x = "Frequencies [cm-1]", y = "Intensities [Debye^2/(amu Å^2)]", data = res)
    plt.xlabel('Wavenumber (cm$^{-1}$)')
    plt.ylabel('Intensity')
    plt.show()
    write_csv_from_dict(res)
    
