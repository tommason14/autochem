__all__ = ['calculate_interaction_energies', 'apply_boltzmann_weightings']

import pandas as pd
from dfply import *
import numpy as np
from ..core.utils import responsive_table

def all_values_are_not_na(column):
    """
    Checks if all values of column are NA, then reverses the boolean.
    So if all values are NA, returns False
    """
    return not all(x for x in np.isnan(column.unique()))

def calculate_interaction_energies(csv, ionic_present=False, software='gamess', pretty_print=False, output=None):
    """
    Calculate interaction energies for ionic clusters from a csv file created
    with this python script; the function assumes that the first subdirectory of
    each path describes a new system, and that purely ionic systems have the
    word 'ionic' in their path, while fragments have the word 'frag' in their
    path.
    Pass in --with-ionic to indicate that a calculation is included that
    includes all ions of the cluster, with neutral/undesired molecules removed.
    """

    @make_symbolic
    def if_else(bools, val_if_true, val_if_false):
        return np.where(bools, val_if_true, val_if_false)
    
    @make_symbolic
    def is_nan(series):
        return np.isnan(series)

    gamess_df = None
    psi4_df = None
    df = pd.read_csv(csv)
    if software == 'gamess':
        gamess_df = df
    else:
        psi4_df = df
                
    if gamess_df is not None:
        if ionic_present:
            data = (gamess_df >>
                mutate(Config = X.Path.str.split('/').str[0]) >>
                mutate(Type = if_else(X.Path.str.contains('frag'), 'frag', 
                                    if_else(X.Path.str.contains('ionic'), 'ionic',
                                            'complex'))) >>
                mutate(HF = X['HF/DFT']) >>
                mutate(SRS = if_else(is_nan(X['MP2/SRS']), X.HF + (1.64 * X.MP2_opp), X['MP2/SRS'])) >>
                mutate(Corr = X.SRS - X.HF) >>
                gather('energy', 'values', [X.HF, X.Corr]) >>
                mutate(energy_type = X.energy + '-' + X.Type) >>
                spread('energy_type', 'values') >>
                group_by(X.Config) >>
                # summing the complexes removes NaN, and need to sum frags anyway for int energy
                summarise(corr_complex = X['Corr-complex'].sum(),
                        corr_ionic = X['Corr-ionic'].sum(),
                        corr_frags = X['Corr-frag'].sum(),
                        hf_complex = X['HF-complex'].sum(),
                        hf_ionic = X['HF-ionic'].sum(),
                        hf_frags = X['HF-frag'].sum()) >>
                mutate(hf_int_kj = (X.hf_complex - X.hf_ionic - X.hf_frags) * 2625.5,
                    corr_int_kj = (X.corr_complex - X.corr_ionic - X.corr_frags) * 2625.5) >>
                mutate(total_int_kj = X.hf_int_kj + X.corr_int_kj)        
            )
            if output is not None:
                data.to_csv(output, index=False)
            if pretty_print:
                data = data.to_dict(orient='list')
                responsive_table(data, strings = [1], min_width=16)
            else:
                print(data)
        else:
            data = (gamess_df >>
                mutate(Config = X.Path.str.split('/').str[0]) >>
                mutate(Type = if_else(X.Path.str.contains('frag'), 'frag', 'complex')) >>
                mutate(HF = X['HF/DFT']) >>
                mutate(SRS = if_else(is_nan(X['MP2/SRS']), X.HF + (1.64 * X.MP2_opp), X['MP2/SRS'])) >>
                mutate(Corr = X.SRS - X.HF) >>
                gather('energy', 'values', [X.HF, X.Corr]) >>
                mutate(energy_type = X.energy + '-' + X.Type) >>
                spread('energy_type', 'values') >>
                group_by(X.Config) >>
                # summing the complexes removes NaN, and need to sum frags anyway for int energy
                summarise(corr_complex = X['Corr-complex'].sum(),
                        corr_frags = X['Corr-frag'].sum(),
                        hf_complex = X['HF-complex'].sum(),
                        hf_frags = X['HF-frag'].sum()) >>
                mutate(hf_int_kj = (X.hf_complex - X.hf_frags) * 2625.5,
                    corr_int_kj = (X.corr_complex - X.corr_frags) * 2625.5)>>
                mutate(total_int_kj = X.hf_int_kj + X.corr_int_kj)
            )
            if output is not None:
                data.to_csv(output, index=False)
            if pretty_print:
                data = data.to_dict(orient='list')
                responsive_table(data, strings = [1], min_width=16)
            else:
                print(data)

    if psi4_df is not None:
        if ionic_present:
            data = (psi4_df >>
                mutate(Config = X.Path.str.split('/').str[0]) >>
                mutate(Type = if_else(X.Path.str.contains('frag'), 'frag', 'complex')) >>
                mutate(HF = X['HF/DFT']) >>
                mutate(SRS = X['HF/DFT'] + 1.64 * X.MP2_opp) >>
                mutate(Corr= X.SRS - X.HF) >>
                gather('energy', 'values', [X.HF, X.Corr]) >>
                mutate(energy_type = X.energy + '-' + X.Type) >>
                spread('energy_type', 'values') >>
                group_by(X.Config) >>
                # summing the complexes removes NaN, and need to sum frags anyway
                # for int energy
                summarise(corr_complex = X['Corr-complex'].sum(),
                        corr_ionic = X['Corr-ionic'].sum(),
                        corr_frags = X['Corr-frag'].sum(),
                        hf_complex = X['HF-complex'].sum(),
                        hf_ionic = X['HF-ionic'].sum(),
                        hf_frags = X['HF-frag'].sum()) >>
                mutate(hf_int_kj = (X.hf_complex - X.hf_ionic - X.hf_frags) * 2625.5,
                    corr_int_kj = (X.corr_complex - X.corr_ionic - X.corr_frags) * 2625.5) >>
                mutate(total_int_kj = X.hf_int_kj + X.corr_int_kj)        
            )
            if output is not None:
                data.to_csv(output, index=False)
            if pretty_print:
                data = data.to_dict(orient='list')
                responsive_table(data, strings = [1], min_width=16)
            else:
                print(data)
        else:
            data = (psi4_df >>
                mutate(Config = X.Path.str.split('/').str[0]) >>
                mutate(Type = if_else(X.Path.str.contains('frag'), 'frag', 'complex')) >>
                mutate(HF = X['HF/DFT']) >>
                mutate(SRS = X['HF/DFT'] + 1.64 * X.MP2_opp) >>
                mutate(Corr= X.SRS - X.HF) >>
                gather('energy', 'values', [X.HF, X.Corr]) >>
                mutate(energy_type = X.energy + '-' + X.Type) >>
                spread('energy_type', 'values') >>
                group_by(X.Config) >>
                # summing the complexes removes NaN, and need to sum frags anyway for int energy
                summarise(corr_complex = X['Corr-complex'].sum(),
                        corr_frags = X['Corr-frag'].sum(),
                        hf_complex = X['HF-complex'].sum(),
                        hf_frags = X['HF-frag'].sum()) >>
                mutate(hf_int_kj = (X.hf_complex - X.hf_frags) * 2625.5,
                    corr_int_kj = (X.corr_complex - X.corr_frags) * 2625.5) >>
                mutate(total_int_kj = X.hf_int_kj + X.corr_int_kj)     
            )
            if output is not None:
                data.to_csv(output, index=False)
            if pretty_print:
                data = data.to_dict(orient='list')
                responsive_table(data, strings = [1], min_width=16)
            else:
                print(data)


def apply_boltzmann_weightings(csv, grouping, output):
    """
    Take in a csv produced from `calculate_interaction_energies` and weight configurations
    according a boltzmann distribution of total energy.
    """
    def bp(series, as_percent = False):
        """
        Takes in energies in Hartrees, produces
        probabilities according to a Boltzmann distribution.
        Also works with group by objects.
        """
        R = 8.3145
        T = 298.15
        h_to_kJ = 2625.5
        series = series * h_to_kJ
        diffs = series - series.min()
        exponent = np.exp((-1 * diffs * 1000) / (R * T))
        summed = exponent.sum()
        if as_percent:
            return (exponent / summed) * 100
        return exponent / summed


    df = pd.read_csv(csv)
    df['complex_total_energy'] = df['hf_complex'] + df['corr_complex']
    df['Groups'] = eval(grouping)
    df['weightings'] = df.groupby('Groups')['complex_total_energy'].apply(bp)
    df['hf_weighted'] = df['hf_int_kj'] * df['weightings']
    df['corr_weighted'] = df['corr_int_kj'] * df['weightings']
    weighted = df.groupby('Groups').agg({'hf_weighted': sum, 'corr_weighted': sum})
    weighted.columns = ['Electrostatics', 'Dispersion']
    print(weighted)
    weighted.reset_index().to_csv(output, index=False)
