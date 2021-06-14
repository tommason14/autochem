import pandas as pd
import numpy as np
from ..core.utils import responsive_table

__all__ = ["calculate_interaction_energies", "apply_boltzmann_weightings"]
HART_TO_KJ = 2625.5

def _calc_srs(row):
    """
    Scaling factors from J. Chem. Phys. 146, 064108 (2017); https://doi.org/10.1063/1.4975326
    """
    _c_os = {
        "cc-pVDZ": 1.752,
        "cc-pVTZ": 1.640,
        "cc-pVQZ": 1.689,
        "aug-cc-pVDZ": 1.372,
        "aug-cc-pVTZ": 1.443,
        "aug-cc-pVQZ": 1.591,
    }
    if row["Basis"] not in _c_os:
        raise AttributeError(
            "Basis set not supported. Code currently scales MP2 interactions for: "
            f"{', '.join(list(_c_os.keys())[:-1])} and {list(_c_os.keys())[-1]}, "
            f"not {row['Basis']}"
        )
    if np.isnan(row["MP2/SRS"]):
        return row["HF"] + (_c_os[row["Basis"]] * row["MP2_opp"])
    return row["MP2/SRS"]

def _bp(series, as_percent=False):
    """
    Takes in energies in Hartrees, produces
    probabilities according to a Boltzmann distribution.
    Also works with group by objects.
    """
    R = 8.3145
    T = 298.15
    series = series * HART_TO_KJ
    diffs = series - series.min()
    exponent = np.exp((-1 * diffs * 1000) / (R * T))
    summed = exponent.sum()
    if as_percent:
        return (exponent / summed) * 100
    return exponent / summed

def _confidence(series):
    """
    95% confidence intervals defined as:
        1.96 * standard deviation from the mean / sqrt(number of items)

    1.96 assumes a normal distribution.

    https://www.itl.nist.gov/div898/handbook/prc/section1/prc14.htm
    http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals_print.html
    """
    return 1.96 * series.std() * (series.size ** -0.5)

def calculate_interaction_energies(
    csv, with_ionic=False, pretty_print=False, output=None
):
    """
    Calculate interaction energies for ionic clusters from a csv file created
    with this python script; the function assumes that the first subdirectory of
    each path describes a new system, and that purely ionic systems have the
    word 'ionic' in their path, while fragments have the word 'frag' in their
    path.
    Pass in --with-ionic to indicate that a calculation is included that
    includes all ions of the cluster, with neutral/undesired molecules removed.
    """

    df = pd.read_csv(csv)
    if with_ionic:
        df["Type"] = np.where(
            df["Path"].str.contains("frag"),
            "frag",
            np.where(df["Path"].str.contains("ionic"), "ionic", "complex"),
        )
        df = df.rename(columns={"HF/DFT": "HF"})
        df["SRS"] = df.apply(_calc_srs, axis=1)
        df["Corr"] = df["SRS"] - df["HF"]
        retain = [x for x in df.columns if x not in ("HF", "Corr")]
        df = df.melt(
            id_vars=retain,
            value_vars=["HF", "Corr"],
            value_name="values",
            var_name="energy",
        )
        df["energy_type"] = df["energy"] + "-" + df["Type"]
        # Create groupby variable (path of each cluster) after pivoting
        df = df.pivot(
            index="Path", columns="energy_type", values="values"
        ).reset_index()
        # Use the index in a new column - ./c2mim_ac/ionic -> c2mim_ac as the config
        df["Config"] = df["Path"].str.split("/").str[1]
        # sum to remove NA values
        df = df.groupby("Config").agg(
            {
                "Corr-complex": "sum",
                "Corr-ionic": "sum",
                "Corr-frag": "sum",
                "HF-complex": "sum",
                "HF-ionic": "sum",
                "HF-frag": "sum",
            }
        )
        # groupby just makes the Config column and index by which rows are collected,
        # reset the index to get the column back
        df = df.reset_index()
        df["HF_int_KJ_mol"] = (
            df["HF-complex"] - df["HF-ionic"] - df["HF-frag"]
        ) * HART_TO_KJ
        df["Corr_int_KJ_mol"] = (
            df["Corr-complex"] - df["Corr-ionic"] - df["Corr-frag"]
        ) * HART_TO_KJ
        df["Total_int_KJ_mol"] = df["HF_int_KJ_mol"] + df["Corr_int_KJ_mol"]
    else:
        df["Type"] = np.where(df["Path"].str.contains("frag"), "frag", "complex")
        df = df.rename(columns={"HF/DFT": "HF"})
        df["SRS"] = df.apply(_calc_srs, axis=1)
        df["Corr"] = df["SRS"] - df["HF"]
        retain = [x for x in df.columns if x not in ("HF", "Corr")]
        df = df.melt(
            id_vars=retain,
            value_vars=["HF", "Corr"],
            value_name="values",
            var_name="energy",
        )
        df["energy_type"] = df["energy"] + "-" + df["Type"]
        # Create groupby variable (path of each cluster) after pivoting
        df = df.pivot(
            index="Path", columns="energy_type", values="values"
        ).reset_index()
        # Use the index in a new column - ./c2mim_ac/ionic -> c2mim_ac as the config
        df["Config"] = df["Path"].str.split("/").str[1]
        # sum to remove NA values
        df = df.groupby("Config").agg(
            {
                "Corr-complex": "sum",
                "Corr-frag": "sum",
                "HF-complex": "sum",
                "HF-frag": "sum",
            }
        )
        # groupby just makes the Config column and index by which rows are collected,
        # reset the index to get the column back
        df = df.reset_index()
        df["HF_int_KJ_mol"] = (df["HF-complex"] - df["HF-frag"]) * HART_TO_KJ
        df["Corr_int_KJ_mol"] = (df["Corr-complex"] - df["Corr-frag"]) * HART_TO_KJ
        df["Total_int_KJ_mol"] = df["HF_int_KJ_mol"] + df["Corr_int_KJ_mol"]
    if pretty_print:
        df_ = df[
            ["Config", "HF_int_KJ_mol", "Corr_int_KJ_mol", "Total_int_KJ_mol"]
        ].to_dict(orient="list")
        responsive_table(df_, strings=[1], min_width=16)
    else:
        print(df)
    if output is not None:
        df.to_csv(output, index=False)


def apply_boltzmann_weightings(csv, grouping, pretty_print=False, output=None):
    """
    Take in a csv produced from `calculate_interaction_energies` and weight configurations
    according a boltzmann distribution of total energy.
    """
    df = pd.read_csv(csv)
    df["Total-complex"] = df["HF-complex"] + df["Corr-complex"]
    df["Groups"] = eval(grouping)
    df["Weightings"] = df.groupby("Groups")["Total-complex"].transform(_bp)
    df["HF-weighted"] = df["HF_int_KJ_mol"] * df["Weightings"]
    df["Corr-weighted"] = df["Corr_int_KJ_mol"] * df["Weightings"]
    weighted_energies = df.groupby("Groups").agg(
        {"HF-weighted": "sum", "Corr-weighted": "sum"}
    )
    weighted_energies = weighted_energies.rename(
        columns={"HF-weighted": "HF", "Corr-weighted": "Corr"}
    ).reset_index()

    hf_confidence_intervals = (
        df.groupby("Groups")["HF-weighted"]
        .apply(_confidence)
        .reset_index(name="HF_confidence_interval")
    )
    corr_confidence_intervals = (
        df.groupby("Groups")["Corr-weighted"]
        .apply(_confidence)
        .reset_index(name="Corr_confidence_interval")
    )
    overall = weighted_energies.merge(hf_confidence_intervals, on="Groups").merge(
        corr_confidence_intervals, on="Groups"
    )
    if pretty_print:
        overall_ = overall.to_dict(orient="list")
        responsive_table(overall_, strings=[1], min_width=16)
    if output is not None:
        overall.to_csv(output, index=False)
