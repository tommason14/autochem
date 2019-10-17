__all__ = ['geodesics',
           'get_h_bonds',
           'get_results_class',
           'parse_results',
           'print_freqs',
           'results_table',
           'search_for_coords',
           'thermochemistry']

from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.thermo import (thermo_data,
                           freq_data_gamess,
                           freq_data_gauss)
from ..core.utils import (read_file,
                          get_files,
                          write_csv_from_dict,
                          write_csv_from_nested,
                          responsive_table)
from ..interfaces.gamess_results import GamessResults
from ..interfaces.orca_results import OrcaResults
from ..interfaces.psi_results import PsiResults
from ..interfaces.gaussian_results import GaussianResults
from .make_files_meta import make_files_from_meta
import os
import re

"""
File: grep_results.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Searches all sub dirs for results
"""

def search_for_coords(dir):
    """
    Recursively searched log/out files of optimisations for a successful
    equilibration- then writes to `equil.xyz`. If unsuccesful, writes to
    `rerun/rerun.xyz`, whilst also creating the corresponding input and
    job file.
    """
    for log in get_files(dir, ('.log', '.out')):
        r = get_results_class(log)
        if r is not None and r.is_optimisation():
            print(f'Searching {r.log}',
                  end=" ")
            r.get_equil_coords()
            print()

def get_results_class(log):
    """Return an instance of the desired class- |GamessResults|, |PsiResults|"""
    log_type = get_type(log)
    logs = {'gamess': GamessResults(log),
            'orca': OrcaResults(log),
            'psi': PsiResults(log),
            'gaussian': GaussianResults(log)}
    return logs.get(log_type, None)


def get_type(filepath):
    """
    Read in file, determine calculation type
    """
    for line in read_file(filepath):
        if 'GAMESS' in line:
            return 'gamess'
        elif 'Psi4' in line or 'PSI4' in line:
            return 'psi'
        elif 'Gaussian' in line:
            return 'gaussian'
        elif 'O   R   C   A' in line:
            return 'orca'


def parse_results(dir, filepath_includes):
    """
    Used internally to parse log files for energies
    """
    output = []
    for log in get_files(dir, ('.out', '.log'), filepath_includes=filepath_includes):
        calc = get_results_class(log)
        filetype = get_type(log)
        try:
            if calc.completed():  # add provision for energies of opts only if equilibrium found
                if not calc.is_hessian():
                    print(log)
                    data = calc.get_data()
                    output.append({'data': data, 'type': filetype})
        except AttributeError:  # if log/out files are not logs of calculations
            continue
    return output


def results_table(dir, file_name, string_to_find):
    """
    Prints energies of all log/out files in current and any sub directories to the screen, with the option of saving to csv.
    """
    # lists are faster to fill than dict values
    # order: files, paths, basis, hf, mp2, mp2_opp, mp2_same
    data = [[], [], [], [], [], [], []]

    output = parse_results(dir, filepath_includes=string_to_find)

    def add_data(data, vals):
        """
        NB: vals has to be a fixed order:
            files, paths, basis, hf, mp2, mp2_opp, mp2_same
        """
        for i, val in enumerate(vals):
            data[i].append(val)
        return data
    
    def remove_column_if_all_na(data):
        return {k: v for k, v in data.items() if not all(val is 'NA' for val in v)}

    for result in output:
        data = add_data(data, result['data'])

    keys = ('File', 'Path', 'Basis', 'HF/DFT', 'MP2/SRS', 'MP2_opp', 'MP2_same')
    
    table_data = {}
    for key, val in zip(keys, data):
        table_data[key] = val 

    table_data = remove_column_if_all_na(table_data)

    responsive_table(table_data, strings=[1, 2, 3], min_width=12)
    write_csv_from_dict(table_data, filename=file_name)


def thermochemistry(dir, string_to_find, mult, temp):
    """
    Returns thermochemical data for all the relevant hessian log files in the given directory and
    subdirectories. Saves to csv file.
    """
    collected = \
        {
            'File': [],
            'ZPVE': [],
            'TC': [],
            'S elec': [],
            'S trans': [],
            'S rot': [],
            'S vib': [],
            'S tot': [],
            'TC - TS': []
        }
    print('Print csv for more info')
    for log in get_files(dir, ('.log', '.out'), filepath_includes=string_to_find):
        r = get_results_class(log)
        try:
            if r.completed():
                if r.is_hessian():
                    res = thermo_data(r.log, mult, temp)
                    res['File'] = r.log
                    for k, v in res.items():
                        collected[k].append(v)
        except AttributeError:
            continue
        except UnicodeDecodeError:
            print(f'{log}- UnicodeDecodeError')
            continue

    # add units to dict keys

    kj = ('ZPVE', 'TC', 'TC - TS')
    jmol = ('S elec', 'S trans', 'S rot', 'S vib', 'S tot')
    kv = list(collected.items())
    collected.clear()
    for k, v in kv:
        if k in kj:
            collected[k + ' [kJ/mol]'] = v
        elif k in jmol:
            collected[k + ' [J/(mol K)]'] = v
        else:
            collected[k] = v
    responsive_table({k:v for k,v in collected.items() if
                      k in ['File','S tot [J/(mol K)]']}, 
                      strings=[1], min_width=10)
    name = write_csv_from_dict(collected, filename='thermo.csv')

def print_freqs(dir):
    """
    Prints frequencies and intensities for jobs in the directory
    passed in. Works for GAMESS/Gaussian jobs.
    """
    for f in os.listdir(dir):
        if f.endswith('log') or f.endswith('out'):
            calc = get_results_class(f)
            if calc.is_hessian():
                if get_type(f) == 'gamess':
                    freq_data_gamess(f, called_by_thermo_code=False) 
                else:
                    freq_data_gauss(f, called_by_thermo_code=False) 
            


def get_h_bonds(dir):
    """
    Searches the current directory for xyz files, then attempts to split them
    into fragments, reporting any intermolecular bonding involving
    hydrogen-bonding atoms less than 2 Å apart, of the correct type and within
    45° of linear. Prints to screen, with the option of saving to csv. If user says yes, writes to hbonds.csv
    """
    print('\n', ' ' * 15, 'HYDROGEN BOND DATA\n')
    output = []
    for file in get_files(dir, ['xyz']):
        path, f = os.path.split(file)
        print('Checking', file[2:])
        mol = Molecule(using=file)
        mol.separate()
        res = mol.find_h_bonds()
        for i in res:
            i.insert(0, f)
            i.insert(1, path)
        output += res
    print()
    if len(output) > 0:
        data = {}
        keys = ('File', 'Path', 'Molecule1', 'Atom1',
                'Molecule2', 'Atom2', 'Length', 'Angle')
        for index, value in enumerate(keys):
            data[value] = [val[index] for val in output]
        responsive_table(data, strings=[1, 2, 3, 4, 5, 6], min_width=6)
        print()
        write_csv_from_nested(output, col_names=(
'File', 'Path', 'Molecule1', 'Atom1', 'Molecule2', 
'Atom2', 'Length', 'Angle'), filename='hbonds.csv')


def file_is_gamess(file):
    """ Check first line of file for 'rungms' string """
    with open(file, 'r') as f:
        return 'rungms' in f.readline()

def geodesics(dir, output):
    """
    Recursively pulls geodesic charges from GAMESS calculations.
    Writes to `charges.csv` if desired
    """

    atom_regex = '^\s[A-Za-z]{1,2}\s*[0-9]*.[0-9]*(\s*-?[0-9]*.[0-9]*){3}$'
    charge_regex = '^\s[A-Za-z]{1,2}(\s*-?[0-9]*.[0-9]*){2}$'

    results = []

    files = get_files(dir, ['log'])
    for logfile in files:
        if file_is_gamess(logfile):
            print(logfile)
            path, filename = os.path.split(logfile)
            inpfile = logfile[:-3] + 'inp'

            res = []
            assigned = []

            for line in read_file(inpfile):
                if re.search(atom_regex, line):
                    sym, atnum, x, y, z = line.split()
                    x, y, z = map(float, (x, y, z))
                    res.append([path, Atom(sym, coords = (x, y, z))]) # new key for each coord
            found = False
            counter = 0
            for line in read_file(logfile):
                if 'NET CHARGES:' in line:
                    found = True
                if 'RMS DEVIATION' in line:
                    break
                if found: 
                    if re.search(charge_regex, line):
                        res[counter].append(float(line.split()[1]))
                        counter += 1
            coordinates = [atom[1] for atom in res]
            mol = Molecule(atoms = coordinates)
            mol.separate()
            for atom, r in zip(mol.coords, res):
                path, _, geodesic_charge = r
                results.append([path, atom.index, atom.symbol, geodesic_charge, 
                atom.x, atom.y, atom.z, f"{mol.fragments[atom.mol]['name']}_{atom.mol}"])

    # nested list (one level) to dict
    data = {}
    keys = ('Path', 'Index', 'Element', 'Geodesic', 'Rx', 'Ry', 'Rz', 'Fragment')
    for index, value in enumerate(keys):
        data[value] = [val[index] for val in results]
    responsive_table(data, strings=[1, 3, 8], min_width=10)

    write_csv_from_nested(results, col_names = keys, filename=output)
