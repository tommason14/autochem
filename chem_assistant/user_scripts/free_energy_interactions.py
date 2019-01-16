__all__ = ['calculate_free_energy_interactions']

# Free energy of interaction between:
# - complex and the constituent ions (pure electrostatics)
# - neutral species and the ionic network
# - can be extended to work with the lithium/sodium

#  (from interaction energies found from single points)
#                     |
# deltaH = Dispersion [kJ/(mol IP)] + [TC(complex) - TC(ionic_network) - TC(neutrals)] * 1/ num_ip 

# deltaS = ([Stot(complex) - Stot(ionic) - Stot(neutrals)] / num_ip) * 298.15 / 1000
# (unit conversion J/mol.K -> kJ/mol)

# Also give the electrostatic free energy (between complex and all ions)
# Use total SRS interaction energies, TC(complex) - sum(TC(ions))
# Stot(complex) - sum(Stot(ions))

# give ratio of free energies between water and cluster to the free energy between cluster and ions

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

def find_e_int(path, csvfile):
    """Gets the dispersion component of the interaction energy (per ion pair) for each configuration, each key of the groups dictionary. This value, when temperature corrected, is the enthalpy of interaction."""
    with open(csvfile, 'r') as f:
        for line in f.readlines()[1:]:
            splitup = line.split(',')
            if splitup[0] == path: #filepath is the first column of csv
                if len(splitup) == 15 or len(splitup) == 16:
                    disp_contribution = float(splitup[6])
                    elec = float(splitup[4]) # neutral species included
                else:
                    disp_contribution = 0.0
                    elec = float(splitup[8]) # just ionic clusters- check index- total mp2 is 
    return disp_contribution, elec
def get_tc_s_tot(d):
    """
    
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
            tc = job[2]
            s_tot = job[-3]
            tc, s_tot = map(float, (tc, s_tot))

        if 'ionic' in job:
            if 'frag' in job: # happens if ionic is ran in the frags dir
                tc = job[2]
                s_tot = job[-4]
                tc, s_tot = map(float, (tc, s_tot))
            else:
                tc = job[2]
                s_tot = job[-3]
                tc, s_tot = map(float, (tc, s_tot))

        if 'frag' in job and 'ionic' not in job:
            if 'neutral' in job:
                tc = job[2]
                s_tot = job[-4]
                tc, s_tot = map(float, (tc, s_tot))
            else:
                tc = job[2]
                s_tot = job[-3]
                tc, s_tot = map(float, (tc, s_tot))
        return tc, s_tot

    int_energy_csvfile = check_user_input('Filename of csv containing interaction energies- created by script', lambda item: item.endswith('.csv'), "Please print a name ending in '.csv'")

    results_dict = {}
    for k, v in d.items():
        sum_frags_s_tot = 0.0
        sum_frags_tc = 0.0
        complex_tc = 0.0
        complex_s_tot = 0.0
        ionic_tc = 0.0
        ionic_s_tot = 0.0
        sum_neutral_tc = 0.0
        sum_neutral_s_tot = 0.0
        for job in v:
            file = job[0] # filepath- ..../hess/frags/water_4/
            if 'hess' in file:
                tc, s_tot = get_results_per_job(job) # now determine type of job
                if 'frags' in file:
                    sum_frags_tc += tc
                    sum_frags_s_tot += s_tot
                    for mol in Molecule.Neutrals:
                        if mol in file:
                            sum_neutral_tc += tc
                            sum_neutral_s_tot += s_tot
                if 'complex' in file:
                    complex_tc = tc
                    complex_s_tot = s_tot
                if 'ionic' in file:
                    ionic_tc = tc
                    ionic_s_tot = s_tot

        # calculate num of ion pairs
        num_ions = 0
        for job in v:
            if 'frag' in job and 'neutral' not in job:
                num_ions += 1
        num_ip = num_ions // 2 # floor division, 5 // 2 = 2


        # dispersion interaction per configuration
        disp, elec = find_e_int(k, int_energy_csvfile)
        disp_int_per_ip = disp / num_ip
        elec_int_per_ip = elec / num_ip

        T = 298.15
        J_TO_KJ = 1000

        dH_neutral = disp_int_per_ip + ((complex_tc - ionic_tc - sum_neutral_tc) / num_ip)
        TdS_neutral = ((complex_s_tot - ionic_s_tot - sum_neutral_s_tot) / num_ip) * T / J_TO_KJ
        dG_neutral = dH_neutral - TdS_neutral
        print(k, dG_neutral) 
  
        dH_elec = elec_int_per_ip + ((complex_tc - sum_frags_tc) / num_ip)
        TdS_elec = ((complex_s_tot - sum_frags_s_tot) / num_ip) * T / J_TO_KJ
        dG_elec = dH_elec - TdS_elec
        print(k, dG_elec) 

        # boltzmann weight it- use the same functions as in other file


    #         results_dict[k] =  {'elec_hf': elec_hf, 'elec_mp2': elec_mp2, 'disp_hf': disp_hf, 'disp_mp2': disp_mp2, 'total_hf': total_hf, 'total_mp2': total_mp2, 'total_mp2_per_ip': total_mp2_per_ip,
    #         'dispersion': dispersion, 'electrostatics': electrostatics} # all the neutral stuff
    #     else:
    #         results_dict[k] = {'total_hf': total_hf, 'total_mp2': total_mp2, 'total_mp2_per_ip': total_mp2_per_ip,'dispersion': dispersion, 'electrostatics': electrostatics}        
                
    # return results_dict, num_ip

def calculate_free_energy_interactions(csv):
    groups = group_files(csv, header = True)
    get_tc_s_tot(groups)