__all__ = ['calculate_free_energy_interactions']

import csv
import re
import math
from ..core.molecule import Molecule
from ..core.utils import check_user_input
# pandas is slow, maybe try something different?

def group_files(csv, header = True):
    """
    Parses a csv file produced by python script
    """

    def split_path(path):
        """
        Returns two strings- one upto the molecule directory, the other further into it.
        Example: c4mim/ac/4/p2/spec/frags/water_4/ --> 
        c4mim/ac/4/p2, spec/frags/water_4

        - Data pre-processing
        """
        upto = 0
        path_split = path.split('/')
        for ind, part in enumerate(path_split):
            if part in ('opt', 'spec', 'hess'):
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
            path, zpve, tc, s_elec, s_trans, s_rot, s_vib, s_tot, tc_ts = line.split(',')
            # split path
            molecule, file = split_path(path)     
            if molecule not in groups:
                groups[molecule] = [[file, zpve, tc, s_elec, s_trans, s_rot, s_vib, s_tot, tc_ts]]
            else:
                groups[molecule].append([file, zpve, tc, s_elec, s_trans, s_rot, s_vib, s_tot, tc_ts])# need to be lists, as later, add on a frag term
    return groups

def calculate_energies(d):
    """
    Calculates interaction energies from single point energy calculations
    INT = E_COMPLEX - SUM(E_ION)
    INT_NEUTRAL SPECIES = E_COMPLEX - E_IONIC_CLUSTER - SUM(E_NEUTRAL_SPECIES)
    """
    def get_results_per_job(job):
        """Returns HF and MP2 energies"""
        file, zpve, tc, s_elec, s_trans, s_rot, s_vib, s_tot, tc_ts = job

        # find type of system
        name = '/'.join(file.split('/')[:-1])
        if 'frags' in name:
            job.append('frag')
            for mol in Molecule.Neutrals:
                if mol in name:
                    job.append('neutral')
        if 'ionic' in name:
            job.append('ionic')
        if 'complex' in name:
            job.append('complex')
        
        # CHANGE FROM HERE
        if 'complex' in job:
            _, _, hf, mp2, _ = job
            hf, mp2 = map(float, (hf, mp2))

        if 'ionic' in job:
            if 'frag' in job: # happens if ionic is ran in the frags dir
                _, _, hf, mp2, _, _ = job
                hf, mp2 = map(float, (hf, mp2))
            else:
                _, _, hf, mp2, _ = job
                hf, mp2 = map(float, (hf, mp2))

        if 'frag' in job and 'ionic' not in job:
            if 'neutral' in job:
                _, _, hf, mp2, _, _ = job
                hf, mp2 = map(float, (hf, mp2))
            else:
                _, _, hf, mp2, _, = job
                hf, mp2 = map(float, (hf, mp2))
        return hf, mp2

    HARTREE_TO_kJ = 2625.5

    purely_ionic = True
    results_dict = {}
    for k, v in d.items():
        sum_frags_hf = 0.0
        sum_frags_mp2 = 0.0
        complex_hf = 0.0
        complex_mp2 = 0.0
        ionic_hf = 0.0
        ionic_mp2 = 0.0
        sum_neutral_hf = 0.0
        sum_neutral_mp2 = 0.0
        for job in v:
            file = job[0]
            if 'spec' in file:
                hf, mp2 = get_results_per_job(job) # now determine type of job
                if 'frags' in file:
                    sum_frags_hf += hf
                    sum_frags_mp2 += mp2
                    for mol in Molecule.Neutrals:
                        if mol in file:
                            sum_neutral_hf += hf
                            sum_neutral_mp2 += mp2
                if 'complex' in file:
                    complex_hf = hf
                    complex_mp2 = mp2
                if 'ionic' in file:
                    ionic_hf = hf
                    ionic_mp2 = mp2


        # calculate num of ion pairs
        num_ions = 0
        for job in v:
            if 'frag' in job and 'neutral' not in job:
                num_ions += 1
        num_ip = num_ions // 2 # floor division, 5 // 2 = 2

        elec_hf = (complex_hf - sum_frags_hf) * HARTREE_TO_kJ 
        elec_mp2 = (complex_mp2 - sum_frags_mp2) * HARTREE_TO_kJ
        disp_hf = 0.0
        disp_mp2 = 0.0
        if ionic_hf != 0.0:
            disp_hf = (complex_hf - ionic_hf - sum_neutral_hf) * HARTREE_TO_kJ
            disp_mp2 = (complex_mp2 - ionic_mp2 - sum_neutral_mp2) * HARTREE_TO_kJ
        total_hf = elec_hf + disp_hf
        total_mp2 = elec_mp2 + disp_mp2
        electrostatics = (total_hf / total_mp2) * 100
        dispersion = total_mp2 - total_hf

        total_mp2_per_ip = total_mp2 / num_ip

        if ionic_hf != 0.0:
            purely_ionic = False
            results_dict[k] =  {'elec_hf': elec_hf, 'elec_mp2': elec_mp2, 'disp_hf': disp_hf, 'disp_mp2': disp_mp2, 'total_hf': total_hf, 'total_mp2': total_mp2, 'total_mp2_per_ip': total_mp2_per_ip,
            'dispersion': dispersion, 'electrostatics': electrostatics} # all the neutral stuff
        else:
            results_dict[k] = {'total_hf': total_hf, 'total_mp2': total_mp2, 'total_mp2_per_ip': total_mp2_per_ip,'dispersion': dispersion, 'electrostatics': electrostatics}        
                
    return results_dict, purely_ionic, num_ip

def calculate_free_energy_interactions(csv):
    groups = group_files(csv, header = True)
    return groups