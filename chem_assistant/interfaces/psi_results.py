from ..core.results import Results

import re

__all__ = ["PsiResults"]


class PsiResults(Results):
    """Class defining the results of a PSI4 calculation."""

    def __init__(self, log):
        super().__init__(log)

    def completed(self):
        complete = False
        for line in self.read():
            if "exiting successfully" in line:
                complete = True
        return complete

    def get_runtype(self):
        """
        Returns runtype. For example, for MP2 single points, the line `energy('mp2')` is used. 
        This method returns the string 'energy'.
        """
        for line in self.read():
            # need regex for energy('mp2') or optimize('scf', dertype='hess') (any number of k-v pairs)
            if re.search("[A-z]*\('[A-z0-9]*'(.?\s*[A-z]*='[A-z]*')*\)", line):
                if re.search("[A-z]*\('[A-z0-9]*'\)", line):  # energy('mp2')
                    return line.split("(")[0]
                else:  # optimize('scf', dertype='hess'......)
                    return line.split("(")[0]  
                    # add to this later, using the collect additional data

    @property
    def method(self):
        """
        Returns energy type. For example, for MP2 single points, the line `energy('mp2')` is used. 
        This method returns the string 'mp2'.
        """
        for line in self.read():
            # need regex for energy('mp2') or optimize('scf', dertype='hess') (any number of k-v pairs)
            if re.search("[A-z]*\('[A-z0-9]*'(.?\s*[A-z]*='[A-z]*')*\)", line):
                if re.search("[A-z]*\('[A-z0-9]*'\)", line):  # energy('mp2')
                    return re.search("[A-z]*\('([A-z0-9]*)'\)", line).group(1)
                # else: #optimize('scf', dertype='hess'......)
                # return line.split('(')[0] #add to this later,

    def is_optimisation(self):
        return self.get_runtype() == "optimize"

    def is_spec(self):
        return self.get_runtype() == "energy"

    def is_hessian(self):
        return self.get_runtype() == "frequency"

    @property
    def multiplicity(self):
        for line in self.read():
            if "Geometry (in Angstrom)" in line:
                return int(line.split()[-1].replace(":", ""))

    def _neutral_homo_lumo(self):
        """
        Finds HOMO-LUMO gap for jobs of singlet multiplicity
        """
        found_region = False
        energies = []
        for line in self.read():
            if "Orbital Energies" in line:
                found_region = True
            if "Final Occupation" in line:
                found_region = False
            if found_region and line.strip() is not "":
                energies.append(line)
        for index, line in enumerate(energies):
            if "Virtual" in line:
                homo_lumo = energies[index - 1 : index + 2]

        homo = float(homo_lumo[0].split()[-1])
        lumo = float(homo_lumo[-1].split()[1])
        return homo, lumo

    def _reduced_homo_lumo(self):
        """
        Finds SOMO-LUMO gap for jobs of doublet multiplicity
        """
        found_singly_occupied = False
        found_virtual = False
        singly = []
        virtual = []
        for line in self.read():
            if "Singly Occupied" in line:
                found_singly_occupied = True
            if "Virtual" in line:
                found_singly_occupied = False
                found_virtual = True
            if "Final Occupation" in line:
                found_virtual = False
            if found_singly_occupied and line.strip() is not "":
                singly.append(line)
            if found_virtual and line.strip() is not "":
                virtual.append(line)
                # save time
                if len(virtual) > 3:
                    break
        somo = float(singly[-1].split()[-1])
        lumo = float(virtual[1].split()[1])
        return somo, lumo

    def _homo_lumo_gap(self):
        hartrees_to_eV = 27.21
        if self.multiplicity == 1:
            homo, lumo = self._neutral_homo_lumo()
        else:
            homo, lumo = self._reduced_homo_lumo()
        homo, lumo = map(lambda x: x * hartrees_to_eV, (homo, lumo))
        gap = lumo - homo
        return homo, lumo, gap

    @property
    def homo_lumo_info(self):
        """
        Prints the HOMO-LUMO gap. Finds SOMO-LUMO if multiplicity is 2.
        Returns `self.multiplicity`, SOMO/HOMO (eV), LUMO (eV) and the gap (eV).
        """

        if self.multiplicity == 1:
            homo, lumo, gap = self._homo_lumo_gap()
            transition = "HOMO-LUMO"
        elif self.multiplicity == 2:
            homo, lumo = self._homo_lumo_gap()  # here homo is somo
            transition = "SOMO-LUMO"
        else:
            print(
                f"Error: Only singlet/doublet multiplicities have been accounted for. Ignoring {self.log}"
            )

        return {
            "File": self.file,
            "Path": self.path,
            "Multiplicity": self.multiplicity,
            "Transition": transition,
            "HOMO/SOMO (eV)": homo,
            "LUMO (eV)": lumo,
            "Gap (eV)": gap,
        }

    def get_data(self):
        """
        Returns job data: filename, filepath, basis set, HF/DFT energy, and MP2 opposite
        and same spin parameters if relevant.
        """
        if self.method == "scf":
            return self._scf_data()

        elif self.method == "mp2":
            return self._mp2_data()

    @property
    def basis(self):
        """
        Returns basis set.
        """
        for line in self.read():
            if re.search("basis\s\w*(\-?\w*){1,2}$", line):
                return line.split()[-1]

    @property
    def total_energy(self):
        """
        Returns total energy, printed for scf calculations.
        """
        total = ""
        for line in self.read():
            if "Total Energy =" in line:
                total = float(line.split("=")[1].strip())
        return total

    def _scf_data(self):
        """
        Return data for scf calculations. 
        Note the NAs returned are because of no MP2 data.
        """
        return (
            self.file,
            self.path,
            self.method,
            self.basis,
            self.total_energy,
            "NA",
            "NA",
            "NA",
        )

    @property
    def hf_energy_for_mp2(self):
        """
        Returns 'reference energy' from MP2 calculations.
        """
        HF = ""
        for line in self.read():
            if "Reference Energy          =" in line:
                HF = float(line.split("=")[1].split()[0].strip())
        return HF

    @property
    def mp2_opp(self):
        """
        Returns MP2 opposite spin energy.
        """
        opp = ""
        for line in self.read():
            if "Opposite-Spin Energy      =" in line:
                opp = float(line.split("=")[1].split()[0].strip())
        return opp

    @property
    def mp2_same(self):
        """
        Returns MP2 same spin energy.
        """
        same = ""
        for line in self.read():
            if "Same-Spin Energy          =" in line:
                same = float(line.split("=")[1].split()[0].strip())
        return same

    def _mp2_data(self):
        """
        Returns data for MP2 calculations: filename, filepath, 
        basis set, hf energy, opp spin energy, same spin energy.
        No MP2 data is returned, but should be calculated instead from the 
        HF and MP2 correlation energies by the user, as coefficients of 
        each spin component will vary depending on the basis set.
        """

        return (
            self.file,
            self.path,
            self.method,
            self.basis,
            self.hf_energy_for_mp2,
            "NA",
            self.mp2_opp,
            self.mp2_same,
        )
