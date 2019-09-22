from .atom import Atom
from .utils import read_file, write_csv_from_dict
import os
import re
import subprocess
import sys

__all__ = ['thermo_data']
# __all__ = ['thermo_data', 'make_ir_spectra']

def get_filetype(file):
    """
    Returns either 'gamess' or 'gauss' depending on filetype of log file
    """
    for line in read_file(file):
        if 'gamess' in line.lower():
            return 'gamess'
        if 'gaussian' in line.lower():
            return 'gauss'   

def thermo_initial_geom_gamess(file):
    """Parses GAMESS hessian calculation log file for the initial geometry"""
    atoms = []
    regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){4}"
    inp = file[:-3] + 'inp'
    with open(inp, "r") as f:
        for line in f.readlines():
            if re.search(regex, line):
                sym, _, x, y, z = line.split()
                x, y, z = map(float, (x, y, z))
                atoms.append(Atom(symbol = sym, coords = (x, y, z)))
    with open('geom.input', 'w') as new:
        for atom in atoms:
            new.write(f"{atom.symbol:5s} {str(atom.atnum):3s} {atom.x:>15.10f} {atom.y:>15.10f} {atom.z:>15.10f} \n")

def thermo_initial_geom_gauss(file):
    """Parses Gaussian frequency calculation log file for the initial geometry. Note that
    coordinates here are stored in .job files by default. If not found, then looks for .inp.
    Only works with xyz coordinates, not z-matrices!!!!"""
    atoms = []
    regex = "[A-Za-z]{1,2}(\s*\D?[0-9]{1,3}\.[0-9]{1,10}){3}"
    inp = file[:-3] + 'job'
    if inp not in os.listdir('.'):
        inp = file[:-3] + 'inp'
    with open(inp, "r") as f:
        for line in f.readlines():
            if re.search(regex, line):
                sym, x, y, z = line.split()
                x, y, z = map(float, (x, y, z))
                atoms.append(Atom(symbol = sym, coords = (x, y, z)))
    with open('geom.input', 'w') as new:
        for atom in atoms:
            new.write(f"{atom.symbol:5s} {str(atom.atnum):3s} {atom.x:>15.10f} {atom.y:>15.10f} {atom.z:>15.10f} \n")

def freq_data_gamess(file, write_freqs_to_file = False):
    """Parses GAMESS hessian calculation log file for the frequency data"""
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
    """Parses Gaussian hessian calculation log file for the frequency data"""
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



def run(file, filetype):
    """
    Calls either thermo-gamess.exe or thermo-gauss.exe depending on which filetype is passed in
    """
    if filetype not in ('gamess', 'gauss'):
        sys.exit("Pass in 'gamess' or 'gauss' to `run`")

    thermo_exe = os.path.join(os.path.dirname(os.path.realpath(__file__)), f'thermo-{filetype}.exe')
    command=f"{thermo_exe}" # file mult temp geom
    if filetype == 'gauss':
        multiplicity=int(input('multiplicity: '))
        temp=float(input('temp (K): '))
    p = subprocess.Popen(command, shell=True, 
    stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True) 
    newline = os.linesep
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

def thermo_data(file):
    """Uses a fortran script to produce thermochemical data for GAMESS Hessian calculations and
    GAUSSIAN frequency calculations- the results produced in the log file have been shown to be 
    inaccurate. At < 300 cm⁻¹, rigid rotor fails, and the fortran code implements hindered rotor."""
    filetype=get_filetype(file)
    if filetype == 'gamess': 
        thermo_initial_geom_gamess(file)
        freq_data_gamess(file, write_freqs_to_file = True)
        run(file, 'gamess')
    if filetype == 'gauss':
        thermo_initial_geom_gauss(file)
        freq_data_gauss(file, write_freqs_to_file = True)
        run(file, 'gauss')
    fort = read_fort()    
    data = grep_data(fort)
    cleanup()
    return data

# def make_ir_spectra(file):
#     """Plot of wavenumber against intensities for vibrations found by diagonalisation of a computed
# hessian matrix"""
#
#     import warnings
#     warnings.filterwarnings("ignore")
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#
#     res = freq_data(file)
#     # now add gaussians... interesting project
#     sns.set_style('darkgrid')
#     sns.lineplot(x = "Frequencies [cm-1]", y = "Intensities [Debye^2/(amu Å^2)]", data = res)
#     plt.xlabel('Wavenumber (cm$^{-1}$)')
#     plt.ylabel('Intensity')
#     plt.show()
#     write_csv_from_dict(res, filename = 'freq.data')
#     
