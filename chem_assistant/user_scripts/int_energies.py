import csv
import re
from chem_assistant import Molecule
# pandas is slow, maybe try something different?

def group_files(csv, header = True):
    """
    Parses a csv file produced by python script
    """
    groups = {}   
    with open(csv, "r") as f:
        if header:
            file_obj = f.readlines()[1:]
        else:
            file_obj = f.readlines

        for line in file_obj:
            file, path, basis, hf, mp2 = line.split(',')
            if path not in groups:
                groups[path] = [[file, basis, hf, mp2]]
            else:
                groups[path].append([file, basis, hf, mp2])# need to be lists, as later, add on a frag term
    return groups

def calculate_energies(d):
    """
    INT = E_COMPLEX - SUM(E_ION)
    INT_NEUTRAL SPECIES = E_COMPLEX - E_IONIC_CLUSTER - SUM(E_NEUTRAL_SPECIES)
    """
    new_dict = {}
    for k, v in d.items():
        if k.endswith('spec'):
            complex_hf = 0.0
            complex_mp2 = 0.0
            sum_frag_hf = 0.0
            sum_frag_mp2 = 0.0
            ionic_hf = 0.0
            ionic_mp2 = 0.0
            sum_neu_hf = 0.0
            sum_neu_mp2 = 0.0
            
            for result in v:
                file, basis, hf, mp2 = result
                name = file.split('_')[0]
                if name in Molecule.Cations or name in Molecule.Anions or name in Molecule.Neutrals or name in ('cation', 'anion', 'neutral'):
                    result.append('frag')
                    if name in Molecule.Neutrals:
                        result.append('neutral')
                # add option for 2IP here- now good for bigger systems
                if 'ionic' in file.lower() or '2ip' in file.lower():
                    result.append('ionic')
            # do the maths
                if 'frag' in result:
                    *_, hf, mp2, _ = result
                    hf, mp2 = map(float, (hf, mp2))
                    sum_frag_hf += hf
                    sum_frag_mp2 += mp2
                # if 2ip in result:
                if 'ionic' in result:
                    *_, hf, mp2, _ = result
                    hf, mp2 = map(float, (hf, mp2))
                    ionic_hf += hf
                    ionic_mp2 += mp2
                if 'neutral' in result:
                    *_, hf, mp2, _ = result
                    hf, mp2 = map(float, (hf, mp2))
                    sum_neu_hf += hf
                    sum_neu_mp2 += mp2
                if 'frag' not in result: # and 2IP not in result
                    # one file per path
                    _, _, hf, mp2 = result 
                    hf, mp2 = map(float, (hf, mp2))
                    complex_hf = hf
                    complex_mp2 = mp2
            elec_hf = (complex_hf - sum_frag_hf) * 2625.5
            elec_mp2 = (complex_mp2 - sum_frag_mp2) * 2625.5

            disp_hf = 0.0
            disp_mp2 = 0.0
            if ionic_hf != 0.0: #then reassign
                disp_hf = (complex_mp2 - ionic_hf - sum_neu_hf) * 2625.5
                disp_mp2 = (complex_mp2 - ionic_mp2 - sum_neu_mp2) * 2625.5

            total_hf = elec_hf + disp_hf # always apply- if no neutral molecules, then disp_hf = 0
            total_mp2 = elec_mp2 + disp_mp2

            electro = (total_hf / total_mp2) * 100
            disp = total_mp2 - total_hf

            if ionic_hf != 0.0:
                new_dict[k] = {} # all the neutral stuff
            else:
                new_dict[k] = {'total_hf': total_hf, 'total_mp2': total_mp2, 'disp': disp, 'electro': electro}

            # separate out purely ionic interaction energies from mixed
    return new_dict

    # new_dict[path] = [[[mol1], [mol2], [mol3]]], int_hf, int_mp2, electro, disp]


def write_csv(data, filename):
    with open(filename, "w") as new:
        writer = csv.writer(new)
        writer.writerow(('Int_MP2 = SRS interaction energies if possible',))
        writer.writerow(('Path', 'Elec Int_HF [kJ/mol]', 'Elec Int_MP2 [kJ/mol]', 'Dispersion [kJ/mol]', '% Electrostatics'))
        for k, v in data.items():
            int_hf, int_mp2, electro, disp = v
            writer.writerow((k, int_hf, int_mp2, electro, disp))
        #could separate ILs with a new line??
        # if cation-anion


def sort_data(data):
    """ 
    Adds a -1 to the first configuration of each IL, if not already there, then sorts the data into alphanumerical order
    """
    collapsed = [[k, v] for k, v in data.items()]
    for index, item in enumerate(collapsed):
        path = item[0]
        *parent, il, spec = path.split('/')
        if not re.match('.*-[0-9]*$', il):
            il = il + '-1'
        parent = '/'.join(parent)
        path = '/'.join((parent, il))
        item[0] = path # needs to be a list of lists to rename the filepath
    sorted_data = sorted(collapsed, key = lambda kv: kv[0])
    sorted_dict  = {}
    for data in sorted_data:
        k, v = data
        sorted_dict[k] = v # as of py 3.6, dicts remain ordered- so no need to implement ordered dict
    return sorted_dict

def assign_molecules(data):
    """ 
    Assign a cation and anion to each path.
    """
    for key in data.keys():
        cation = ''
        anion = ''
        vals = key.split('/')
        for val in vals:
            # different names for the same anion
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
                energies.append((path, data['total_mp2']))
            sorted_vals = sorted(energies, key = lambda tup: tup[1]) #sort on the energies
            for index, val in enumerate(sorted_vals, 1): #start index at 1
                path, energy = val
                groups[cation][anion][path]['rank'] = index
                groups[cation][anion][path]['deltaE'] = groups[cation][anion][path]['total_mp2'] / sorted_vals[0][1] # minimum energy

    for cation in groups:
        for anion in groups[cation]:
            print(cation, anion)
            for path, data in groups[cation][anion].items():
                print(path, data['rank'])

    return groups

def calculate_interaction_energies(csv):
    data = group_files(csv)
    new_data = calculate_energies(data)
    sorted_data = sort_data(new_data)
    assigned = assign_molecules(sorted_data)
    ranked = rank_configs(assigned)
    correct = False
    while not correct:
        filename = input('Filename: ')
        if not filename.endswith('.csv'):
            print("Please end the filename in '.csv'")
        else:
            correct = True
    write_csv(groups, filename)

calculate_interaction_energies('1ip_energies.csv')