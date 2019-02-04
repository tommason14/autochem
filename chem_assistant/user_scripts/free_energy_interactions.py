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
from ..core.utils import check_user_input, sort_data
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
                if 'opt' in cell or 'spec' in cell:
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
    disp_contribution = 0.0
    elec = 0.0
    if path[0] == '.':
        config = path[2:]
    else:
        config = path
    with open(csvfile, 'r') as f:
        for line in f.readlines()[1:]:
            splitup = line.split(',')
            if splitup[0] == config : #filepath is the first column of csv
                if len(splitup) == 15 or len(splitup) == 16:
                    disp_contribution = float(splitup[6])
                    elec = float(splitup[4]) # neutral species included
                else:
                    disp_contribution = 0.0
                    elec = float(splitup[8]) # just ionic clusters- check index- total mp2 is elec
    return disp_contribution, elec


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


def calc_free_energies(d):
    """
    
    """
    int_energy_csvfile = check_user_input('Filename of csv containing interaction energies- created by script', lambda item: item.endswith('.csv'), "Please print a name ending in '.csv'")

    results = {}
    for config, v in d.items():
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
        disp, elec = find_e_int(config, int_energy_csvfile)
        disp_int_per_ip = disp / num_ip
        elec_int_per_ip = elec / num_ip

        T = 298.15
        J_TO_KJ = 1000

        dG_neutral = 0.0

        if sum_neutral_tc != 0.0:
            dH_neutral = disp_int_per_ip + ((complex_tc - ionic_tc - sum_neutral_tc) / num_ip)
            TdS_neutral = ((complex_s_tot - ionic_s_tot - sum_neutral_s_tot) / num_ip) * T / J_TO_KJ
            dG_neutral = dH_neutral - TdS_neutral
  
        dH_elec = elec_int_per_ip + ((complex_tc - sum_frags_tc) / num_ip)
        TdS_elec = ((complex_s_tot - sum_frags_s_tot) / num_ip) * T / J_TO_KJ
        dG_elec = dH_elec - TdS_elec
    
        dG_total = dG_elec + dG_neutral
        
        results[config] = {'dH_neutral': dH_neutral, 'dH_elec': dH_elec,
        'TdS_neutral': TdS_neutral, 'TdS_elec': TdS_elec, 'dG_neutral': dG_neutral, 'dG_elec': dG_elec, 'dG_total': dG_total}
        
    return results


    #         results_dict[k] =  {'elec_hf': elec_hf, 'elec_mp2': elec_mp2, 'disp_hf': disp_hf, 'disp_mp2': disp_mp2, 'total_hf': total_hf, 'total_mp2': total_mp2, 'total_mp2_per_ip': total_mp2_per_ip,
    #         'dispersion': dispersion, 'electrostatics': electrostatics} # all the neutral stuff
    #     else:
    #         results_dict[k] = {'total_hf': total_hf, 'total_mp2': total_mp2, 'total_mp2_per_ip': total_mp2_per_ip,'dispersion': dispersion, 'electrostatics': electrostatics}        
                
    # return results_dict, num_ip

def assign_molecules_from_dict_keys(data):
    """ 
    Assign a cation and anion to each path.
    """
    for key in data.keys():
        cation = ''
        anion = ''
        vals = key.split('/')
        for val in vals:
            # different names for the same anion
            if val == 'ch':
                val = 'choline'
            if val == 'ac':
                val = 'acetate'
            if val == 'h2po4':
                val = 'dhp' # in Molecules.Anions
            if val == 'mesylate':
                val = 'mes'
            if val in Molecule.Cations:
                cation = val
            elif val in Molecule.Anions:
                anion = val
        data[key]['cation'] = cation
        data[key]['anion'] = anion
    return data


def rank_configs(data):
    """
    Ranks each configuration according to its interaction energies.
    """
    KJ_TO_J = 1000
    R = 8.3145
    T = 298.15

    # create groups
    cations = sorted(set([v['cation'] for v in data.values()]))
    anions  = sorted(set([v['anion']  for v in data.values()]))
    # add to groups in alphabetical order

    groups = {}
    for c in cations:
        groups[c] = {}
        for a in anions:
            groups[c][a] = {}

    # add data to correct bin
    for k, v in data.items():
        for cation in cations:
            for anion in anions:
                if cation == v['cation'] and anion == v['anion']:
                    groups[cation][anion][k] = v
    
    for cation in groups:
        for anion in groups[cation]:
            energies = []
            for path, data in groups[cation][anion].items():
                if data['dG_neutral'] != 0.0: 
                    energies.append((path, data['dG_neutral']))
                else:
                    energies.append((path, data['dG_elec']))
            sorted_vals = sorted(energies, key = lambda tup: tup[1]) #sort on the energies
            min_energy = sorted_vals[0][1]
            for index, val in enumerate(sorted_vals, 1): #start index at 1
                path, energy = val
                if data['dG_neutral'] != 0.0:
                    ddG = groups[cation][anion][path]['dG_neutral'] - min_energy
                else:
                    ddG = groups[cation][anion][path]['dG_elec'] - min_energy
                groups[cation][anion][path]['rank'] = index
                groups[cation][anion][path]['ddG'] = ddG
                groups[cation][anion][path]['boltzmann_factor'] =\
                math.exp((-1 * KJ_TO_J * ddG) / (R * T))
                # missing a division somewhere- when we use the BF to weight the average

    # change the order of the paths of each config in the groups dict for each cation-anion pair, by their rank

    ordered_dict = {}
    for cat in groups:
        ordered_dict[cat] = {}
        for an in groups[cat]:
            ordered_dict[cat][an] = {}
            lst = [(k, v) for k, v in groups[cat][an].items()]
            sorted_lst = sorted(lst, key = lambda kv: kv[1]['rank'])
            for kv in sorted_lst:
                k, v = kv
                ordered_dict[cat][an][k] = v 
    return ordered_dict

def write_csv(data, filename):

    neutral_included = False
    for cation in data:
        for anion in data[cation]:
            for config in data[cation][anion]:
                if data[cation][anion][config]['dG_neutral'] != 0.0:
                    neutral_included = True
                    break # check once only

    def calc_boltz_ave(d, neu):
        boltz_ave_neu = 0.0
        boltz_ave_elec = 0.0
        boltz_ave_tot = 0.0

        KJ_TO_J = 1000
        R = 8.3145
        T = 298.15
        lookup = 'dG_neutral' if neu else 'dG_total'
        
        total_prob = 0.0
        for config in d:
            total_prob += math.exp((-1 * d[config][lookup]) / (R * T))
            d[config]['boltzmann_factor'] = d[config]['boltzmann_factor'] / total_prob
        for config in d:
            if neu:
                boltz_ave_neu += d[config]['dG_neutral'] * d[config]['boltzmann_factor']
                boltz_ave_elec += d[config]['dG_elec'] * d[config]['boltzmann_factor']
                boltz_ave_tot += d[config]['dG_total'] * d[config]['boltzmann_factor']
            else:
                boltz_ave_tot += d[config]['dG_total'] * d[config]['boltzmann_factor']
        
        if neu:
            return boltz_ave_elec, boltz_ave_neu, boltz_ave_tot
        else:
            return boltz_ave_tot
            # needs adding only once per cat-an

    if neutral_included:
        col_names = ('Path', 'Cation', 'Anion',
        'ΔH Electrostatics [kJ/(mol IP)]', 'ΔH neutral [kJ/(mol IP)]',
        'TΔS Electrostatics [kJ/(mol IP)]', 'TΔS neutral [kJ/(mol IP)]',
        'ΔG Electrostatics [kJ/(mol IP)]', 'ΔG neutral [kJ/(mol IP)]', 
        'ΔG Total [kJ/(mol IP)]','ΔΔG Neutral [kJ/(mol IP)]', 'Rank', 'Boltzmann Weighting', 
        'BW ΔG Electrostatics [kJ/(mol IP)]', 'BW ΔG Neutral [kJ/(mol IP)]', 'BW ΔG Total [kJ/(mol IP)]')

        variables = ('dH_elec', 'dH_neutral', 'TdS_elec', 'TdS_neutral', 'dG_elec', 'dG_neutral', 'dG_total', 'ddG', 'rank', 'boltzmann_factor') 
    else:
        col_names = ('Path', 'Cation', 'Anion', 
        'ΔH Total [kJ/(mol IP)]','TΔS Total [kJ/(mol IP)]',
        'ΔG Total [kJ/(mol IP)]','ΔΔG Total [kJ/(mol IP)]', 'Rank', 
        'Boltzmann Weighting', 'BW ΔG Total [kJ/(mol IP)]')

        variables = ('dH_elec', 'TdS_elec', 'dG_total', 'ddG', 'rank', 'boltzmann_factor')
    
    def update_csv(path, cat, an, d, variables):
        locals().update(d) #create variables
        lst = [path, cat, an]
        for key in variables:
            lst.append(locals()[key])
        return lst

    with open(filename, "w", encoding = 'utf-8-sig') as new:
        writer = csv.writer(new)
        # writer.writerow(('Int_MP2 = SRS interaction energies if possible',))
        writer.writerow(col_names)
        for cation in data:
            for anion in data[cation]:
                if neutral_included:
                    boltz_ave_elec, boltz_ave_neu, boltz_ave_tot =\
                    calc_boltz_ave(data[cation][anion], neutral_included)
                else:
                    boltz_ave_tot =\
                    calc_boltz_ave(data[cation][anion], neutral_included)

                for index, value in enumerate(data[cation][anion].items()):
                    path, d = value
                    numbers = update_csv(path, cation, anion, d, variables)
                    if index == 0:
                        if neutral_included:
                            numbers = numbers + [boltz_ave_elec, boltz_ave_neu, boltz_ave_tot]
                        else:
                            numbers = numbers + [boltz_ave_tot]
                    writer.writerow(numbers)
 


def calculate_free_energy_interactions(csv):
    groups = group_files(csv, header = True)
    res = calc_free_energies(groups)
    sorted_data = sort_data(res)
    assigned = assign_molecules_from_dict_keys(sorted_data)
    ranked = rank_configs(assigned)
    filename = check_user_input('Filename of output', lambda item: item.endswith('.csv'), "Please enter a filename ending in '.csv'")

    write_csv(ranked, filename)
