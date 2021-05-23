from .atom import Atom
from .utils import read_file, write_csv_from_dict
import os
from glob import glob
import re
import subprocess
import sys

__all__ = ["thermo_data", "freq_data_gamess", "freq_data_gauss"]


def get_filetype(file):
    """
    Returns either 'gamess' or 'gauss' depending on filetype of log file
    """
    for line in read_file(file):
        if "gamess" in line.lower():
            return "gamess"
        if "gaussian" in line.lower():
            return "gauss"


def write_geom_input(atoms):
    """
    Writes 'geom.input' from the list of |Atom| objects passed in.
    """
    with open("geom.input", "w") as new:
        for atom in atoms:
            new.write(
                f"{atom.symbol:5s} {int(atom.atomic_number):3} {atom.x:>15.10f} {atom.y:>15.10f} {atom.z:>15.10f} \n"
            )


def write_freq_out_file(results):
    """
    Writes freq.out using the frequencies of the `results` dictionary.
    Note all of the list is written to the file, so any removal of 
    rotations and translations must occur before using this function.
    """
    with open("freq.out", "w") as output:
        for i in results["Frequencies [cm-1]"]:
            output.write(f"{i:.3f}\n")


def thermo_initial_geom_gamess(file):
    """
    Parses GAMESS output for the initial geometry.
    Takes the nuclear coordinates from the 'coord 0 vib 0' run,
    and converts from Bohrs to angstroms.
    """
    atoms = []
    BOHR_TO_ANG = 0.529177
    regex = "\s+[0-9]+\s+[A-z]+(\s*-?[0-9]+\.[0-9]+){3}"
    found = False
    for line in read_file(file):
        if re.search("VIB\s+0", line):
            found = True
        if found and "CPU" in line:
            break
        if found and re.search(regex, line):
            _, sym, x, y, z = line.split()
            x, y, z = map(lambda v: float(v) * BOHR_TO_ANG, (x, y, z))
            atoms.append(Atom(symbol=sym, coords=(x, y, z)))
    write_geom_input(atoms)


def rm_additional_rots_and_trans(results):
    """
    Removes the first 6 rotations and translations, to leave
    the required 3N vibrations.
    """
    for key, value in results.items():
        results[key] = value[6:]
    return results


def thermo_initial_geom_gauss(file):
    """
    Parses Gaussian frequency calculation log file for the initial 
    geometry. Note that coordinates here are stored in .job files by     
    default. Only works with xyz coordinates, not z-matrices.
    """
    atoms = []
    regex = "\s+[A-z]{1,2}(\s+-?[0-9]+\.[0-9]+){3}"
    # found_freq = False
    found_coords = False
    for line in read_file(file):
        if "Symbolic Z-matrix:" in line:
            found_coords = True
        if line == "\n":
            found_coords = False
        if found_coords and re.search(regex, line):
            sym, x, y, z = line.split()
            x, y, z = map(float, (x, y, z))
            atoms.append(Atom(symbol=sym, coords=(x, y, z)))

    write_geom_input(atoms)


def freq_data_gamess(file):
    """Parses GAMESS hessian log files for frequency data"""
    regex = (
        "[0-9]{1,9}?\s*[0-9]{1,9}\.[0-9]{1,9}\s*[A-Za-z](\s*[0-9]{1,9}\.[0-9]{1,9}){2}$"
    )
    found_region = False
    modes = []
    freqs = []
    ints = []
    with open(file, "r") as f:
        for line in f:
            if "MODE FREQ(CM**-1)  SYMMETRY  RED. MASS  IR INTENS." in line:
                found_region = True
            if line == "\n":
                found_region = False
            if found_region:
                if re.search(regex, line):
                    mode, vib, *_, intensity = line.split()
                    mode = int(mode)
                    vib, intensity = map(float, (vib, intensity))
                    modes.append(mode)
                    freqs.append(vib)
                    ints.append(intensity)

    results = {
        "Modes": modes,
        "Frequencies [cm-1]": freqs,
        "Intensities [Debye^2/(amu Å^2)]": ints,
    }  # keys used as headers for csv

    results = rm_additional_rots_and_trans(results)
    write_freq_out_file(results)
    return results


def freq_data_gauss(file):
    """Parses Gaussian frequency log files for frequency data"""
    freqs = []
    ints = []
    with open(file, "r") as f:
        for line in f:
            if "Frequencies --" in line:
                freqs += line.split()[2:]
            if "IR Inten    --" in line:
                ints += line.split()[3:]

    freqs = [float(i) for i in freqs]
    ints = [float(i) for i in ints]
    modes = [i for i in range(len(freqs))]
    results = {
        "Modes": modes,
        "Frequencies [cm-1]": freqs,
        "Intensities [Debye^2/(amu Å^2)]": ints,
    }  # keys used as headers for csv

    results = rm_additional_rots_and_trans(results)
    write_freq_out_file(results)

    return results


def run(file, mult, temp):
    """
    Calls thermo.exe with geom.input and freq.out written to the same 
    directory.
    """
    thermo_exe = os.path.join(os.path.dirname(os.path.realpath(__file__)), "thermo.exe")
    p = subprocess.Popen(
        thermo_exe,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    newline = os.linesep
    commands = ["y", "y", "y", mult, temp]
    p.communicate(newline.join(commands))


def read_fort():
    with open("fort.10", "r") as f:
        fort = [line for line in f.readlines()]
    return fort


def grep_data(fort):
    """Collects desired data from fort.10"""
    data = {}
    lookup = {
        "TC h": "TC",
        "S elec": "S elec",
        "S trans": "S trans",
        "S rot": "S rot",
        "S vib": "S vib",
        "Stot": "S tot",
        "TC-": "TC - TS",
    }
    for line in fort:
        if re.search("ZPVE.*=.*kJ", line):
            data["ZPVE"] = line.split()[2]
        for k, v in lookup.items():
            if k in line:
                data[v] = line.split()[-1]
    return data


def cleanup():
    os.system("rm fort.10 moments geom.input freq.out")


def setup_and_run_fortran_script(file, mult, temp):
    """
    Runs fortran script to produce 'fort.10' files etc...
    """
    filetype = get_filetype(file)
    if filetype == "gamess":
        thermo_initial_geom_gamess(file)
        freq_data_gamess(file)
        run(file, mult, temp)
    if filetype == "gauss":
        thermo_initial_geom_gauss(file)
        freq_data_gauss(file)
        run(file, mult, temp)


def thermo_data(file, mult, temp):
    """
    Uses a fortran script to produce thermochemical data for GAMESS 
    Hessian calculations and GAUSSIAN frequency calculations- the results 
    produced in the GAMESS files have been shown to be 
    inaccurate. At < 300 cm⁻¹, rigid rotor fails, and the fortran code 
    implements hindered rotor.
    """
    setup_and_run_fortran_script(file, mult, temp)
    fort = read_fort()
    data = grep_data(fort)
    cleanup()
    return data
