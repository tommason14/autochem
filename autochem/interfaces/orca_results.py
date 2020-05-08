from ..core.utils import read_file, write_xyz
from ..core.results import Results
from ..core.periodic_table import PeriodicTable as PT
from ..core.atom import Atom

import re
import os
import subprocess

__all__ = ["OrcaResults"]


class OrcaResults(Results):
    """
    Class for obtaining results from Orca simulations. This class     
    requires a log file to be read.
    """

    def __init__(self, log):
        super().__init__(log)

    def get_runtype(self):
        """
        Can have multiple runs in one file i.e. opt freq
        Looks for one or both of 'opt' and 'freq'.
        If not found, defaults to single point calculation.
        """
        opt = False
        freq = False
        td = False
        if "opt" in self.user_commands:
            opt = True
        if "freq" in self.user_commands:
            freq = True
        if "td" in self.user_commands:
            td = True
        if opt and not freq and not td:
            return "opt"
        if opt and freq:
            return "opt-freq"
        if td and opt:
            return "es-opt"
        if td:
            return "vert"
        if freq and not opt:
            return "freq"
        return "spec"

    def completed(self):
        for line in self.eof(0.05):
            if "****ORCA TERMINATED NORMALLY****" in line:
                return True
        return False

    def is_optimisation(self):
        return "opt" in self.get_runtype()

    def is_spec(self):
        return self.get_runtype() == "spec"

    def is_hessian(self):
        return "freq" in self.get_runtype()

    @property
    def user_commands(self):
        """
        Returns the ! ... line of the input file.
        """
        for line in self.read():
            if "> !" in line:
                return line.lower()

    @property
    def title(self):
        """
        Returns xyz file with no extension. Used when writing new coords
        """
        for line in self.read():
            if "The coordinates will be read from file" in line:
                return line.split()[-1].rsplit(".")[0]

    @property
    def is_dft(self):
        """
        Used internally to decide if dft energies should be collected.
        """
        for line in self.read():
            if "Density Functional     Method          .... DFT" in line:
                return True
        return False

    @property
    def method(self):
        """
        Returns method used in calculation.
        """
        for line in self.read():
            # dft
            if "Exchange Functional    Exchange" in line:
                return line.split()[-1]
            # HF
            elif "Ab initio Hamiltonian  Method" in line:
                return line.split()[-1].split("(")[0]
            # MP2
            # elif ....

    @property
    def num_atoms(self):
        for line in self.read():
            if "Number of atoms" in line:
                return int(line.split()[-1])

    def get_equil_coords(self):
        coords = []
        found_coords = False
        found_equil = False
        regex = "^\s+[A-z]+(\s+-?[0-9]+\.[0-9]+){3}$"

        for line in self.eof(0.5):
            if "THE OPTIMIZATION HAS CONVERGED" in line:
                found_equil = True
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                found_coords = True
            if line is "\n":
                found_coords = False
            if found_coords and re.search(regex, line):
                if len(coords) == self.num_atoms:
                    coords = []
                sym, x, y, z = line.split()
                coords.append(Atom(sym, coords=[x, y, z]))

        if found_equil:
            print("Found equilibrium!")
            newdir = os.path.join(self.path, "spec")
            if not os.path.isdir(newdir):
                os.mkdir(newdir)
            write_xyz(coords, os.path.join(newdir, f"{self.title}-equil.xyz"))
        else:
            if len(coords) > 0:
                print(
                    "Equilibrium not found. Needs resubmitting.",
                    f"Coords stored in {self.path}/rerun/{self.title}.xyz",
                )
                newdir = os.path.join(self.path, "rerun")
                if not os.path.isdir(newdir):
                    os.mkdir(newdir)
                write_xyz(coords, os.path.join(newdir, f"{self.title}-rerun.xyz"))
            else:
                print("No iterations were cycled through!")

    @property
    def basis(self):
        """
        Returns basis set.
        """
        for line in self.read():
            if "Your calculation utilizes the basis:" in line:
                return line.split()[-1]

    @property
    def total_energy(self):
        """
        Returns total energy, printed for scf calculations.
        """
        for line in self.read():
            if "Total Energy       :" in line:
                return float(line.split()[3].strip())

    @property
    def final_single_point_energy(self):
        """
        Returns the energies printed for single points.
        """
        for line in self.eof(0.2):
            if "FINAL SINGLE POINT ENERGY" in line:
                return float(line.split()[-1])

    @property
    def sp_data(self):
        """
        Returns data for HF/DFT single point calculations.
        Note the NAs returned are because of no MP2 spin parameters.
        """
        return (
            self.file,
            self.path,
            self.method,
            self.basis,
            self.final_single_point_energy,
            "NA",
            "NA",
        )

    def get_data(self):
        """
        Returns job data: filename, filepath, basis set, HF/DFT energy, and MP2 opposite
        and same spin parameters if relevant.
        """
        if self.is_spec:
            return self.sp_data

    @property
    def multiplicity(self):
        """
        Return multiplicity.
        """
        for line in self.read():
            if "Multiplicity" in line:
                return int(line.split()[-1])

    def _homo_lumo(self):
        """
        Finds HOMO/LUMO orbitals, and ORCA prints the
        energies in eV, so no need for conversion.
        """
        homo = ""
        lumo = ""
        regex = r"^\s+[0-9]+(\s+-?[0-9]+.[0-9]+){3}"
        found = False
        for line in self.read():
            if "ORBITAL ENERGIES" in line:
                found = True
            if found:
                if re.search(regex, line):
                    line = line.split()
                    if line[1] != "0.0000":  # occupancy
                        homo = line[-1]
                    else:
                        lumo = line[-1]
                        break
        homo, lumo = map(float, (homo, lumo))
        return homo, lumo

    def _homo_lumo_gap(self):
        homo, lumo = self._homo_lumo()
        gap = lumo - homo
        return homo, lumo, gap

    @property
    def homo_lumo_info(self):
        """
        Prints the HOMO-LUMO gap. Finds SOMO-LUMO if multiplicity is 2.
        Returns `self.multiplicity`, SOMO/HOMO (eV), LUMO (eV) and the gap (eV).
        """

        if self.multiplicity == 1:
            transition = "HOMO-LUMO"
        else:
            transition = "SOMO-LUMO"
        
        homo, lumo, gap = self._homo_lumo_gap()
        
        return {
            "File": self.file,
            "Path": self.path,
            "Multiplicity": self.multiplicity,
            "Transition": transition,
            "HOMO/SOMO (eV)": homo,
            "LUMO (eV)": lumo,
            "Gap (eV)": gap,
        }

    # Vibrational analysis

    @property
    def frequencies(self):
        """
        Orca removes rotations/vibrations before printing.
        """
        vibs = []
        found = False
        regex = r"^\s+[0-9]+:(\s+-?[0-9]+.[0-9]+){2}\s+\((\s+-?[0-9]+.[0-9]+){3}\)"
        for line in self.read():
            if "Mode    freq (cm**-1)" in line:
                found = True
            if line is "\n":
                found = False
            if found and re.search(regex, line):
                vibs.append(float(line.split()[1]))
        return vibs

    @property
    def intensities(self):
        """
        Orca removes rotations/vibrations before printing
        """
        ints = []
        regex = r"^\s+[0-9]+:(\s+-?[0-9]+.[0-9]+){2}\s+\((\s+-?[0-9]+.[0-9]+){3}\)"
        for line in self.read():
            if "Mode    freq (cm**-1)" in line:
                found = True
            if line is "\n":
                found = False
            if found and re.search(regex, line):
                ints.append(float(line.split()[2]))
        return ints

    # TD-DFT Excited states

    @property
    def td_dft_wavelengths(self):
        """
        Returns a nested list of wavelengths per iteration
        """
        waves = []
        waves_per_iter = []
        found = False
        regex = "^\s+[0-9]+(\s+-?[0-9]+\.[0-9]+){7}$"
        for line in self.read():
            if "TRANSITION ELECTRIC" in line:
                found = True
            if found:
                if re.search(regex, line):
                    wave = float(line.split()[2])
                    waves_per_iter.append(wave)
            if line is "\n":
                found = False
                if len(waves_per_iter) > 0:
                    waves.append(waves_per_iter)
                    waves_per_iter = []
        return waves

    @property
    def td_dft_intensities(self):
        """
        Returns a nested list of intensities
        """
        ints = []
        ints_per_iter = []
        found = False
        regex = "^\s+[0-9]+(\s+-?[0-9]+\.[0-9]+){7}$"
        for line in self.read():
            if "TRANSITION ELECTRIC" in line:
                found = True
            if found:
                if re.search(regex, line):
                    intensity = float(line.split()[4])
                    ints_per_iter.append(intensity)
            if line is "\n":
                found = False
                if len(ints_per_iter) > 0:
                    ints.append(ints_per_iter)
                    ints_per_iter = []
        return ints

    @property
    def td_dft_transition_energies(self):
        """
        Returns a nested list of energies of each transition, converted 
        from cm-1 to eV
        """
        inverse_cm_to_ev = 1 / 8065.6
        vals = []
        vals_per_iter = []
        found = False
        regex = "^\s+[0-9]+(\s+-?[0-9]+\.[0-9]+){7}$"
        for line in self.read():
            if "TRANSITION ELECTRIC" in line:
                found = True
            if found:
                if re.search(regex, line):
                    val = float(line.split()[1]) * inverse_cm_to_ev
                    vals_per_iter.append(val)
            if line is "\n":
                found = False
                if len(vals_per_iter) > 0:
                    vals.append(vals_per_iter)
                    vals_per_iter = []
        return vals
