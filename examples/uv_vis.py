#!/usr/bin/env python3

"""
File: uv_vis.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: This pulls UV-Vis data from Orca and Gaussian
log files. Note that this script does not account for fluorescence,
but only because of the NA given for the root. This can easily be 
changed.
"""

from chem_assistant import GaussianResults, OrcaResults, get_log_type
from glob import glob


def results(logfile):
    """
    Returns the results class of Orca and Gaussian log files.
    If log is of a different type, returns None.
    """
    _types = {"gaussian": GaussianResults, "orca": OrcaResults}

    for key, val in _types.items():
        if get_log_type(logfile) is key:
            return val(logfile)


with open("uv_vis.csv", "w") as f:
    f.write(
        "Config,"
        "Root,"
        "Iteration,"
        "Transition Energies (eV),"
        "Wavelength (nm),"
        "Intensity (au)\n"
    )
    for logfile in glob("**/*log", recursive=True):
        log = results(logfile)
        if log is None:
            continue  # if not gaussian/orca job
        for iteration, data in enumerate(
            zip(
                log.td_dft_wavelengths,
                log.td_dft_intensities,
                log.td_dft_transition_energies,
            ),
            1,
        ):
            waves, ints, energies = data
            for wave, intensity, energy in zip(waves, ints, energies):
                f.write(f"{logfile},NA,{iteration},{energy},{wave},{intensity}\n")
