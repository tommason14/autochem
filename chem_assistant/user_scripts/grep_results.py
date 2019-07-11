__all__ = ['create_extra_jobs',
           'get_h_bonds',
           'get_results_class',
           'parse_results',
           'results_table',
           'search_for_coords',
           'thermochemistry']

from ..core.molecule import Molecule
from ..core.thermo import thermo_data
from ..core.utils import (read_file,
                          get_files,
                          write_csv_from_dict,
                          write_csv_from_nested,
                          responsive_table)
from ..interfaces.gamess_results import GamessResults
from ..interfaces.psi_results import PsiResults
from ..interfaces.gaussian_results import GaussianResults
from .make_files_meta import make_files_from_meta
import os

"""
File: grep_results.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Searches all sub dirs for results
"""
"""
- Create opt dir or spec or hess-  when creating the jobs (make another dir)
TODO:
- If geom_opt has completed, extract equilibrium coords:
  save in parent folder as equil.xyz & if no spec directory available,
  create it. If no equil.xyz in spec dir, then copy equil.xyz over to
  spec dir, and make files accordingly.

Get results.
...
...
...
...
Create single points for completed geom_opts? [y/n]
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
        if r is not None:
            if r.completed():
                if r.is_optimisation():
                    print(f'{r.log}:\nFinding equilibrium coordinates...',
                          end=" ")
                    r.get_equil_coords()
                    print()
            else:
                print(f"{log}: Not completed\n")


def create_extra_jobs(dir):
    """
    Using `equil.xyz` files found from a previous search, create
    additional files
    """

    def is_complex():
        complex = False
        while not complex:
            try:
                user_choice = input(
                    ("Are these systems made of ",
                     "more than one molecule? [Y/N] "))
            except ValueError:
                print("Please enter 'Y' or 'N'")
            if user_choice.upper() not in ('Y', 'N'):
                print("Please enter 'Y' or 'N'")
            else:
                complex = True
        return complex

    def add_to_meta(scomp, complexes):
        with open('meta.py', "r") as f:
            lines = [line for line in f]
        lines.insert(-1, f"s.supercomp = '{scomp}'\n")
        if complexes:
            lines[-1] = lines[-1][:-1] + ', is_complex = True)'
        with open('meta.py', "w") as writer:
            for line in lines:
                writer.write(line)

    def copy_meta(file_to_copy, scomp):
        """
        Give the path and filename of meta.py file to copy over to any
        directory containing an `equil.xyz` file
        """
        complex = is_complex()
        parent = os.getcwd()
        for path, _, files in os.walk('.'):
            for file in files:
                if file == 'equil.xyz':
                    os.chdir(path)
                    os.system(f'cp {file_to_copy} .')
                    add_to_meta(scomp, complex)
                    os.chdir(parent)

    def get_supercomp():
        sc = False
        while not sc:
            try:
                sc_choice = int(input(
                    """Which supercomputer do you want to use?

1. Raijin
2. Magnus
3. Massive

Choice: """))
            except ValueError:
                print('Please enter a number between 1 and 3')
            if sc_choice not in range(1, 4):
                print('Please enter a number between 1 and 3')
            else:
                sc = True

        supercomps = {1: 'rjn', 2: 'mgs', 3: 'mas'}
        scomp = supercomps[sc_choice]
        return scomp

    # os.path.dirname(__file__) = user_scripts
    template_dir = os.path.join(os.path.dirname(__file__),
                                os.pardir,
                                'templates/meta_files')
    predefined = {
        1: os.path.join(template_dir, 'gamess/spec/meta.py'),
        2: os.path.join(template_dir, 'gamess/freq/meta.py'),
        3: os.path.join(template_dir, 'psi/spec/meta.py')
    }
    good = False
    while not good:
        try:
            user_choice = int(input(
                """Use a predefined file or one of your own:

1. Predefined
2. Own meta.py

Choice: """))
        except ValueError:
            print('Please choose 1 or 2')
        if user_choice in (1, 2):
            good = True
        else:
            print('Please choose a number, 1 or 2')
    if user_choice == 1:
        done = False
        while not done:
            try:
                choice = int(input(
                    """Which type of file do you want to make?

1. Gamess single point energy
2. Gamess hessian
3. Psi4 single point energy

Choice: """))
            except ValueError:
                print('Please enter a number between 1 and 3')
            if choice not in range(1, 4):
                print('Please enter a number between 1 and 3')
            else:
                done = True

        # add desired meta.py
        meta_file = predefined[choice]
        scomp = get_supercomp()
        copy_meta(meta_file, scomp)
    else:
        user_path = input(("Give path to your own meta.py ",
                           "(relative from the directory you are running ",
                           "this script from)"))
        user_given = os.path.join(os.getcwd(), user_path)
        scomp = get_supercomp()
        copy_meta(user_path, scomp)
        # load script from os.path.join(os.getcwd(), p)
        # if no meta.py there, ask again

    # run meta.py recursively
    make_files_from_meta('.')


def get_results_class(log):
    """Return an instance of the desired class- |GamessResults|, |PsiResults|"""
    log_type = get_type(log)
    logs = {'gamess': GamessResults(log),
            'psi': PsiResults(log),
            'gaussian': GaussianResults(log)}
    return logs.get(log_type, None)


def get_type(filepath):
    """
    Read in file, determine calculation type
    """
    calc = ''
    for line in read_file(filepath):
        if 'GAMESS' in line:
            calc = 'gamess'
            break
        elif 'PSI4' in line:
            calc = 'psi'
            break
        elif 'Gaussian' in line:
            calc = 'gaussian'
    return calc


def parse_results(dir):
    """
    Used internally to parse log files for energies
    """
    output = []
    for log in get_files(dir, ('.out', '.log')):
        calc = get_results_class(log)
        filetype = get_type(log)
        try:
            if calc.completed():  # add provision for energies of opts only if equilibrium found
                if not calc.is_hessian():
                    print(f"Searching through {calc.log}")
                    data = calc.get_data()
                    output.append({'data': data, 'type': filetype})
            # else:
                # print(f'{calc.log}: Incomplete') # don't care about knowing which files are
                # incomplete
        except AttributeError:  # if log/out files are not logs of calculations
            continue
    return output


def results_table(dir):
    """
    Prints energies of all log/out files in current and any sub directories to the screen, with the option of saving to csv.
    """
    output = parse_results(dir)
    table_data = {'File': [],
                  'Path': [],
                  'Basis': [],
                  'HF/DFT': [],
                  'MP2/SRS': [],
                  'MP2_opp': [],
                  'MP2_same': []}

    def add_to_table(table_data, f, p, b, hf, mp2, mp2_opp, mp2_same):
        table_data['File'].append(f)
        table_data['Path'].append(p)
        table_data['Basis'].append(b)
        table_data['HF/DFT'].append(hf)
        table_data['MP2/SRS'].append(mp2)
        table_data['MP2_opp'].append(mp2_opp)
        table_data['MP2_same'].append(mp2_same)
        return table_data
    
    def remove_column_if_all_na(data):
        return {k: v for k, v in data.items() if not all(val is 'NA' for val in v)}

    for result in output:
        if result['type'] == 'psi':
            f, p, b, hf, mp2_opp, mp2_same = result['data']
            mp2 = 'NA'
            table_data = add_to_table(
                table_data, f, p, b, hf, mp2, mp2_opp, mp2_same)
        else:
            f, p, b, hf, mp2 = result['data']
            mp2_opp = 'NA'
            mp2_same = 'NA'
            table_data = add_to_table(
                table_data, f, p, b, hf, mp2, mp2_opp, mp2_same)

    table_data = remove_column_if_all_na(table_data)

    responsive_table(table_data, strings=[1, 2, 3], min_width=12)
    name = write_csv_from_dict(table_data, filename = 'energies.csv')
    # for use in other calculations (chem_assist -r uses this name)
    return name


def thermochemistry(dir):
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
            'S tran': [],
            'S rot': [],
            'S vib': [],
            'S tot': [],
            'TC - TS': []
        }

    for log in get_files(dir, ('.log', '.out')):
        r = get_results_class(log)
        try:
            if r.completed():
                if r.is_hessian():
                    print(f'Thermo data for {r.log}')
                    res = thermo_data(r.log)  # run fortran script
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
    jmol = ('S elec', 'S tran', 'S rot', 'S vib', 'S tot')
    kv = list(collected.items())
    collected.clear()
    for k, v in kv:
        if k in kj:
            collected[k + ' [kJ/mol]'] = v
        elif k in jmol:
            collected[k + ' [J/(mol K)]'] = v
        else:
            collected[k] = v
    name = write_csv_from_dict(
        collected, return_name=True, filename='thermo.csv')
    return name


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
                'Molecule2', 'Atom2', 'Length (Å)', 'Angle (°)')
        for index, value in enumerate(keys):
            data[value] = [val[index] for val in output]
        responsive_table(data, strings=[1, 2, 3, 4, 5, 6], min_width=6)
        print()
        write_csv_from_nested(output, col_names=('File', 'Path', 'Molecule1', 'Atom1',
                                                 'Molecule2', 'Atom2', 'Length (Å)', 'Angle (°)'), filename='hbonds.csv')
