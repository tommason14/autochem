#!/usr/bin/env python3

import csv
import sys

def group_files(csv, header = True):
    """
    Parses a csv file produced by python script
    """

    def split_path(path):
        """
        Returns two strings- one upto the files directory, the other    
        further into it.
        Example: opts/amps/monomer/files/a1/ -->
        opts/amps/monomer/, a1

        - Data pre-processing
        """
        upto = 0
        path_to_mol = ''
        path_to_each_calc = ''
        path_split = path.split('/')
        for ind, part in enumerate(path_split):
            if part == 'files':
                upto = ind
                break
        if upto != 0:
            path_to_mol = path_split[0: upto]
            path_to_each_calc = path_split[upto + 1:]

        path_to = '/'.join(path_to_mol)
        path_after_run = '/'.join(path_to_each_calc)
        return path_to, path_after_run


    groups = {}
    with open(csv, "r") as f:
        if header:
            file_obj = f.readlines()[1:]
        else:
            file_obj = f.readlines

        for line in file_obj:
            # disregard any opts
            for cell in line.split(','):
                if 'opt' in cell:
                    continue
            file, path, basis, dft, _ = line.split(',')
            # split path
            molecule, path_to_file = split_path(path)
            if molecule not in groups:
                groups[molecule] = [[file, path, basis, dft]]
            else:
                groups[molecule].append([file, path, basis, dft])

    ordered = {}
    copy = groups.copy()
    for k, lst in copy.items():
        lst = sorted(lst, key = lambda val: float(val[3]))
        min_hf = float(lst[0][3])
        for i, val in enumerate(lst):
            rank = i + 1
            diff = (float(val[3]) - min_hf) * 2625.5
            val.append(rank)
            val.append(diff)
        ordered[k] = lst
    return ordered


def write_csv_from_dict(data, filename, col_names):
    """Write to file from dictionary"""
    with open(filename, "w", encoding = 'utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(col_names)
        for v in data.values():
            for i in v:
                writer.writerow(i)


def rank_systems(csv):
    """Rank configs according to HF/DFT energy
    """
    data = group_files(csv)
    write_csv_from_dict(data, csv, col_names = ('File', 'Path', 'Basis', 'DFT [Eh]', 'Rank', 'ΔE [kJ/mol]'))

rank_systems(sys.argv[1])
