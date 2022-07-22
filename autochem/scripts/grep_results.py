from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.thermo import thermo_data, freq_data_gamess, freq_data_gauss
from ..core.utils import (
    check_user_input,
    eof,
    get_files,
    list_of_dicts_to_one_level_dict,
    read_file,
    responsive_table,
    write_csv_from_dict,
    write_csv_from_nested,
)
from ..interfaces.gamess_results import GamessResults
from ..interfaces.orca_results import OrcaResults
from ..interfaces.psi_results import PsiResults
from ..interfaces.gaussian_results import GaussianResults
import os
import re
import sys
import pandas as pd

__all__ = [
    "charges",
    "energies",
    "energy_table",
    "file_as_results_class",
    "get_h_bonds",
    "homo_lumo_gaps",
    "nmr_shieldings",
    "print_freqs",
    "print_freqs_to_csv",
    "search_for_coords",
    "thermochemistry",
]


def search_for_coords(dir):
    """
    Recursively searched log/out files of optimisations for a successful
    equilibration- then writes to `equil.xyz`. If unsuccesful, writes to
    `rerun/rerun.xyz`, whilst also creating the corresponding input and
    job file.
    """

    def checked_before(r):
        """
        check that no rerun or spec for that particular file has been created before
        """
        reruns = False
        equils = False
        if "rerun" in os.listdir(r.path):
            reruns = any(f"{r.title}" in f for f in os.listdir(f"{r.path}/rerun"))
        if "spec" in os.listdir(r.path):
            equils = any(f"{r.title}" in f for f in os.listdir(f"{r.path}/spec"))
        return reruns or equils

    for log in get_files(dir, (".log", ".out")):
        r = file_as_results_class(log)
        if r is not None and r.is_optimisation():
            if not checked_before(r):
                print(f"Searching {r.log}")
                r.get_equil_coords()
                print()


def file_as_results_class(log):
    """Return an instance of the desired class- |GamessResults|, |PsiResults|"""
    log_type = get_type(log)
    logs = {
        "gamess": GamessResults(log),
        "orca": OrcaResults(log),
        "psi": PsiResults(log),
        "gaussian": GaussianResults(log),
    }
    return logs.get(log_type, None)


def get_type(filepath):
    """
    Read in file, determine calculation type
    """
    for line in read_file(filepath):
        if "GAMESS" in line:
            return "gamess"
        elif "Psi4" in line or "PSI4" in line:
            return "psi"
        elif "Gaussian" in line:
            return "gaussian"
        elif "O   R   C   A" in line:
            return "orca"


def need_gauss_energy(calc):
    """
    Returns True is there is an energy to be pulled from a Gaussian file.
    Required as Gaussian hessian calculations can be preceeded by optimisations,
    so the calc.is_hessian() call is irrelevant.
    """
    return (
        isinstance(calc, GaussianResults) and calc.is_optimisation() or calc.is_spec()
    )


def energies(dir, filepath_includes):
    """
    Used internally to parse log files for energies
    """
    output = []
    for log in get_files(dir, (".out", ".log"), filepath_includes=filepath_includes):
        calc = file_as_results_class(log)
        filetype = get_type(log)
        try:
            if (
                calc.completed()
            ):  # add provision for energies of opts only if equilibrium found
                if not calc.is_hessian() or need_gauss_energy(calc):
                    print(log)
                    data = calc.get_data()
                    output.append({"data": data, "type": filetype})
        except AttributeError:  # if log/out files are not logs of calculations
            continue
    return output


def energy_table(dir, file_name, string_to_find=None, autosave=None):
    """
    Prints energies of all log/out files in current and any sub directories to the screen,
    with the option of saving to csv.
    """
    # lists are faster to fill than dict values
    # order: file, path, method, basis, hf, mp2, mp2_opp, mp2_same
    data = [[], [], [], [], [], [], [], []]
    # at some point, will make this a dictionary, loads clearer that way.

    output = energies(dir, filepath_includes=string_to_find)

    def add_data(data, vals):
        """
        NB: vals has to be a fixed order:
            files, paths, basis, hf, mp2, mp2_opp, mp2_same
        """
        for i, val in enumerate(vals):
            data[i].append(val)
        return data

    def remove_column_if_all_na(data):
        return {k: v for k, v in data.items() if not all(val == "NA" for val in v)}

    for result in output:
        data = add_data(data, result["data"])

    keys = (
        "File",
        "Path",
        "Method",
        "Basis",
        "HF/DFT",
        "MP2/SRS",
        "MP2_opp",
        "MP2_same",
    )
    table_data = {}
    for key, val in zip(keys, data):
        table_data[key] = val

    table_data = remove_column_if_all_na(table_data)

    if len(table_data) == 0:
        sys.exit("No optimisations or single points found")

    responsive_table(table_data, strings=[1, 2, 3], min_width=12)
    write_csv_from_dict(table_data, filename=file_name, autosave=autosave)


def homo_lumo_gaps(dir, output, string_to_find=None, autosave=None):
    """
    Returns HOMO-LUMO or SOMO-LUMO gaps for each single point calculation
    found in any subdirectory. Currently restricted to single points for
    simplicity, but can probably be extended to optimisations if needed-
    would have to check the log files first.
    """
    info = []
    for log in get_files(dir, (".out", ".log"), filepath_includes=string_to_find):
        calc = file_as_results_class(log)
        filetype = get_type(log)
        try:
            if calc.completed() and calc.is_spec():
                data = calc.homo_lumo_info
                info.append(data)
        except AttributeError:  # if log/out files are not logs of calculations
            continue
    if len(info) == 0:
        sys.exit("Error: No single points found")
    info = list_of_dicts_to_one_level_dict(info)
    responsive_table(info, strings=[1, 2, 4])
    write_csv_from_dict(info, filename=output, autosave=autosave)
    return info


def thermochemistry(dir, string_to_find, mult, temp, output, autosave=None):
    """
    Returns thermochemical data for all the relevant hessian log files in the given directory and
    subdirectories. Saves to csv file.
    """
    collected = {
        "File": [],
        "Method": [],
        "Basis": [],
        "Temperature [K]": [],
        "Multiplicity given": [],
        "ZPVE": [],
        "TC": [],
        "S tot": [],
        "S elec": [],
        "S trans": [],
        "S rot": [],
        "S vib": [],
        "TC - TS": [],
    }
    print("Print csv for more info")
    for log in get_files(dir, (".log", ".out"), filepath_includes=string_to_find):
        r = file_as_results_class(log)
        try:
            if r.completed():
                if r.is_hessian():
                    res = thermo_data(r.log, mult, temp)
                    res["File"] = r.log
                    res["Method"] = r.method
                    res["Basis"] = r.basis
                    res["Temperature [K]"] = temp
                    res["Multiplicity given"] = mult

                    for k, v in res.items():
                        collected[k].append(v)
        except AttributeError:
            continue
        except UnicodeDecodeError:
            print(f"{log}- UnicodeDecodeError")
            continue

    # add units to dict keys

    kj = ("ZPVE", "TC", "TC - TS")
    jmol = ("S tot", "S elec", "S trans", "S rot", "S vib")
    kv = list(collected.items())
    collected.clear()
    for k, v in kv:
        if k in kj:
            collected[k + " [kJ/mol]"] = v
        elif k in jmol:
            collected[k + " [J/(mol K)]"] = v
        else:
            collected[k] = v

    responsive_table(
        {
            k: v
            for k, v in collected.items()
            if k
            in ("File", "Temperature [K]", "Multiplicity given", "S tot [J/(mol K)]")
        },
        strings=[1],
        min_width=10,
    )
    name = write_csv_from_dict(collected, filename=output, autosave=autosave)


def print_freqs(dir, output, string_to_find=None, autosave=None):
    """
    Writes frequencies and intensities of GAMESS/Gaussian frequency calculations
    to a csv. Works recursively through the file system.
    """
    data = {}
    data["File"] = []
    data["Frequencies"] = []
    data["Intensities"] = []
    for file in get_files(dir, ["log", "out"], filepath_includes=string_to_find):
        if "slurm" not in file:
            calc = file_as_results_class(file)
            if calc.is_hessian():
                data["Frequencies"] += calc.frequencies
                data["Intensities"] += calc.intensities
                data["File"] += [calc.log] * len(calc.frequencies)
    responsive_table(data, strings=[1])
    write_csv_from_dict(data, filename=output, autosave=autosave)


def print_freqs_to_csv(dir):
    """
    Writes a comma-separated file of frequencies and intensities for each frequency calculation
    in the current directory
    """
    for f in os.listdir(dir):
        if f.endswith("log") or f.endswith("out"):
            calc = file_as_results_class(f)
            if calc.is_hessian():
                print(f"Extracting freqs from {calc.log}")
                with open(f"{calc.basename}.ir.data.csv", "w") as f:
                    f.write("Frequencies,Intensities\n")
                    for val in zip(calc.frequencies, calc.intensities):
                        freq, intensity = val
                        f.write(f"{freq},{intensity}\n")


def get_h_bonds(dir, output=None, string_to_find=None, autosave=None):
    """
    Searches the current directory for xyz files, then attempts to split them
    into fragments, reporting any intermolecular bonding involving
    hydrogen-bonding atoms less than 2 Å apart, of the correct type and within
    45° of linear. Prints to screen, with the option of saving to csv. If user says yes, writes to hbonds.csv
    """

    def can_cast_as_float(item):
        try:
            item = float(item)
        except ValueError:
            return False
        return True

    distance = check_user_input(
        "Distance (Å) [2]",
        lambda item: can_cast_as_float(item) or item == "",
        "Please enter a number",
    )

    if distance == "":
        distance = 2.0
    else:
        distance = float(distance)

    print("\n", " " * 15, "HYDROGEN BOND DATA\n")
    output = []
    files = [
        f
        for f in get_files(dir, ["xyz"], filepath_includes=string_to_find)
        if f.count("/") == 1
    ]
    for file in files:
        path, f = os.path.split(file)
        print("Checking", file[2:])
        mol = Molecule(using=file)
        mol.separate()
        res = mol.find_h_bonds(distance)
        for i in res:
            i.insert(0, f)
            i.insert(1, path)
        output += res
    print()
    if len(output) > 0:
        data = {}
        keys = (
            "File",
            "Path",
            "Molecule1",
            "Atom1",
            "Molecule2",
            "Atom2",
            "Length",
            "Angle",
        )
        for index, value in enumerate(keys):
            data[value] = [val[index] for val in output]
        responsive_table(data, strings=[1, 2, 3, 4, 5, 6], min_width=9)
        print()
        write_csv_from_nested(
            output,
            col_names=(
                "File",
                "Path",
                "Molecule1",
                "Atom1",
                "Molecule2",
                "Atom2",
                "Length",
                "Angle",
            ),
            filename="hbonds.csv",
            autosave=autosave,
        )


def file_is_gamess(file):
    """Check first line of file for 'rungms' string"""
    with open(file, "r") as f:
        return "rungms" in f.readline()


def file_is_gaussian(file):
    """Check for the word Gaussian in the first 5 lines"""
    for number, line in enumerate(read_file(file)):
        if "gaussian" in line.lower():
            return True
        if number == 5:
            break
    return False


def charges(dir, output, string_to_find=None, autosave=None):
    """
    Recursively pulls geodesic charges from GAMESS calculations.
    Pulls mulliken charges from Gaussian calculations.
    Writes to `charges.csv` if desired
    """
    results = []
    files = get_files(dir, ["log"], filepath_includes=string_to_find)
    for logfile in files:
        if file_is_gaussian(logfile):
            res = []
            print(logfile)
            atom_regex = "^\s?[A-z]{1,2}(\s+-?[0-9]+\.[0-9]+){3}"
            charge_regex = "^\s+[0-9]+\s+[A-z]{1,2}\s+-?[0-9]+\.[0-9]+"
            #     1  C   -0.122119
            for line in read_file(logfile):
                if re.search(atom_regex, line):
                    sym, x, y, z = line.split()
                    x, y, z = map(float, (x, y, z))
                    res.append(
                        [logfile, Atom(sym, coords=(x, y, z))]
                    )  # new key for each coord
            found = False
            counter = 0
            for line in eof(logfile, 0.20):
                if "Mulliken charges:" in line:
                    found = True
                if "Sum of Mulliken charges" in line:
                    break
                if found:
                    if re.search(charge_regex, line):
                        res[counter].append(float(line.split()[-1]))
                        counter += 1
            coordinates = [atom[1] for atom in res]
            mol = Molecule(atoms=coordinates)
            mol.separate()
            for atom, r in zip(mol.coords, res):
                path, _, charge = r
                try:
                    results.append(
                        [
                            path,
                            atom.index,
                            atom.symbol,
                            charge,
                            "NA",  # need to pull out standard devation...
                            atom.x,
                            atom.y,
                            atom.z,
                            f"{mol.fragments[atom.mol]['name']}_{atom.mol}",
                        ]
                    )
                except KeyError:
                    results.append(
                        [
                            path,
                            atom.index,
                            atom.symbol,
                            charge,
                            "NA",  # need to pull out standard devation...
                            atom.x,
                            atom.y,
                            atom.z,
                            "NA",
                        ]
                    )

        if file_is_gamess(logfile):
            atom_regex = "^\s[A-Za-z]{1,2}\s*[0-9]*.[0-9]*(\s*-?[0-9]*.[0-9]*){3}$"
            charge_regex = "^\s[A-Za-z]{1,2}(\s*-?[0-9]*.[0-9]*){2}$"
            print(logfile)
            inpfile = logfile[:-3] + "inp"

            res = []

            for line in read_file(inpfile):
                if re.search(atom_regex, line):
                    sym, atnum, x, y, z = line.split()
                    x, y, z = map(float, (x, y, z))
                    res.append(
                        [logfile, Atom(sym, coords=(x, y, z))]
                    )  # new key for each coord
            found = False
            counter = 0
            for line in read_file(logfile):
                if "NET CHARGES:" in line:
                    found = True
                if "RMS DEVIATION" in line:
                    break
                if found:
                    if re.search(charge_regex, line):
                        _, charge, esd = line.split()
                        charge, esd = map(float, (charge, esd))
                        res[counter].append(charge)
                        res[counter].append(esd)
                        counter += 1
            coordinates = [atom[1] for atom in res]
            mol = Molecule(atoms=coordinates)
            mol.separate()
            for atom, r in zip(mol.coords, res):
                path, _, geodesic_charge, esd = r
                try:
                    results.append(
                        [
                            path,
                            atom.index,
                            atom.symbol,
                            geodesic_charge,
                            esd,
                            atom.x,
                            atom.y,
                            atom.z,
                            f"{mol.fragments[atom.mol]['name']}_{atom.mol}",
                        ]
                    )
                except KeyError:
                    results.append(
                        [
                            path,
                            atom.index,
                            atom.symbol,
                            geodesic_charge,
                            esd,
                            atom.x,
                            atom.y,
                            atom.z,
                            "NA",
                        ]
                    )
    # nested list (one level) to dict
    data = {}
    keys = (
        "Path",
        "Index",
        "Element",
        "Charge",
        "Estimated Stdev",
        "Rx",
        "Ry",
        "Rz",
        "Fragment",
    )
    for index, value in enumerate(keys):
        data[value] = [val[index] for val in results]
    responsive_table(data, strings=[1, 3, 9], min_width=10)

    write_csv_from_nested(results, col_names=keys, filename=output, autosave=autosave)


def nmr_shieldings(dir, output, string_to_find=None, autosave=None):
    shifts = []
    for log in get_files(dir, (".out", ".log"), filepath_includes=string_to_find):
        calc = file_as_results_class(log)
        try:
            if calc.completed() and calc.is_spec():
                data = calc.isotropic_nmr_shielding_constants.assign(
                    Path=calc.path, File=calc.file
                )[["Path", "File", "Index", "Element", "Shielding"]]
                shifts.append(data)
        except AttributeError:  # if log/out files are not logs of calculations
            continue
    if len(shifts) == 0:
        sys.exit("Error: No single points found")
    shifts = pd.concat(shifts, ignore_index=True)
    responsive_table(shifts, strings=[1, 2, 4])
    write_csv_from_dict(shifts.to_dict("list"), filename=output, autosave=autosave)
