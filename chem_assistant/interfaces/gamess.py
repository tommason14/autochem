from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.settings import Settings, read_template, dict_to_settings
from ..core.job import Job
from ..core.periodic_table import PeriodicTable as PT
from ..core.sc import Supercomp
from ..core.utils import consecutive, sort_elements, write_xyz

from os import chdir, mkdir, getcwd, system, walk, listdir
from os.path import exists, join, dirname

__all__ = ["GamessJob"]


class GamessJob(Job):
    # Note the job scripts require the supercomputer to be entered, such as:

    # >>> j = GamessJob(using = 'file.xyz')
    # >>> j.supercomp = 'raijin'
    """Class for creating GAMESS input files and job scripts. 

    The input files generated default to geometry optimisations at the SRS-MP2/cc-pVDZ level of theory. This is easily changed by creating a |Settings| object and adding parameters, using the following syntax.
    >>> s = Settings()
    >>> s.input.contrl.runtyp = 'energy'
    >>> s.input.basis.gbasis = 'CCT'
    >>> s.input.mp2.scsopo = 1.64
    >>> j = GamessJob(using = '../xyz_files/ch_ac.xyz', settings = s)
    This yields the following result:
         $SYSTEM MEMDDI=0 MWORDS=500 $END
         $CONTRL ICHARG=0 ISPHER=1 MAXIT=200 RUNTYP=ENERGY SCFTYP=RHF $END
         $STATPT NSTEP=500 $END
         $SCF DIIS=.TRUE. DIRSCF=.TRUE. FDIFF=.FALSE. $END
         $BASIS GBASIS=CCT $END
         $MP2 CODE=IMS SCSOPO=1.64 SCSPAR=0.0 SCSPT=SCS $END
         $DATA
        ch_ac
        C1
         H 1.0
         C 6.0
         N 7.0
         O 8.0
         $END
         C     6.0   -6.27719   -2.98190   -7.44828
         H     1.0   -7.26422   -2.51909   -7.63925
         H     1.0   -5.81544   -3.32313   -8.39516 
    If FMO (Fragment Molecular Orbital) calculations are desired, pass the keyword argument ``fmo``, set to *True*, along with a settings object with a parameter of ``nfrags``:
        >>> s = Settings()
        >>> s.nfrags = 4
        >>> job = GamessJob(using = 'file.xyz', fmo = True, settings = s)

    The class creates different subdirectories for every molecule in the system.
    Using the class with `frags_in_subdir` set to true produces:
        - a `complex` subdirectory- one per xyz
        - an `ionic` subdirectory for every complex with the neutral species and/or single atom ions removed
            - for water inclusion, N2 inclusion, alkali metal inclusion
        - a `frags` subdirectory for every fragment of the complex

    The names of files created default to the type of calculation: optimisation (opt), single point
energy (spec) or hessian matrix calculation for thermochemical data and vibrational frequencies
(hess). If a different name is desired, pass a string with the ``filename`` parameter, with no extension. The name will be used for both input and job files.
        >>> job = GamessJob(using = 'file.xyz', fmo = True, filename = 'benzene')

    This command produces two files, benzene.inp and benzene.job.

    Files are placed in a subdirectory of their own name. So when creating optimisation files, files are placed in opt:
        .
        └── opt
            ├── opt.inp
            └── opt.job

    To group fragments together, just pass in a `grouped` argument in the settings file:

        >>> sett = Settings()
        >>> sett.grouped = 'water-chloride'

    This will group water and chloride fragments together, to avoid having lots of nodes with a
    small number of atoms assigned to them. 

    """

    def __init__(
        self,
        using=None,
        fmo=False,
        frags_in_subdir=False,
        settings=None,
        filename=None,
        is_complex=False,
        run_dir=None,
        bonds_to_split=None,
    ):
        # Also read in bonds to split from settings object
        self.fmo = fmo  # Boolean
        self.filename = filename
        self.defaults = read_template("gamess.json")  # settings object
        if settings is not None:
            self.user_settings = settings.as_dict()
            self.merged = self.defaults.merge(settings)  # merges inp, job data
            self.input = self.merged.input
            self.job = self.merged.job
        else:
            self.input = self.defaults.input

        # Split each molecule on a bond?
        if bonds_to_split is None:
            if "bonds_to_split" in self.merged:
                bonds_to_split = self.merged.bonds_to_split

        super().__init__(using, user_settings=settings, bonds_to_split=bonds_to_split)

        if "/" in using:
            # say using = ../xyz_files/file.xyz --> file
            self.title = using.split("/")[-1][:-4]
        else:
            self.title = using[:-4]
        self.xyz = using

        self.create_complex_dir_if_required(is_complex, frags_in_subdir)

        self.made_run_dir = False
        if run_dir is not None:
            self.made_run_dir = True

        self.create_inp()
        self.create_job()
        self.place_files_in_dir()

        if frags_in_subdir:
            self.create_inputs_for_fragments(complex_is_fmo=self.fmo)

    def create_complex_dir_if_required(self, is_complex, make_frags):
        self.is_complex = is_complex
        if make_frags and not is_complex:
            self.is_complex = True

    def determine_fragments(self):
        if self.fmo:
            self.mol.separate()
            fmo_data = self.fmo_formatting()
            self.input.fmo = fmo_data
            self.input.fmoprp.maxit = 200
            self.input.gddi.ngroup = len(self.mol.fragments)

    def order_header(self):
        if self.fmo:
            desired = [
                "SYSTEM",
                "CONTRL",
                "GDDI",
                "STATPT",
                "SCF",
                "BASIS",
                "FMO",
                "FMOPRP",
            ]  # mp2/dft after
        else:
            desired = ["SYSTEM", "CONTRL", "STATPT", "SCF", "BASIS"]

        self.header = []
        for i in desired:
            for line in self.unordered_header:
                item = line.split()[0][1:]
                if item == i:
                    self.header.append(line)
        # add in additional commands after, given as sett.input.blah = 'blah'
        else:
            for line in self.unordered_header:
                item = line.split()[0][1:]
                if item not in desired:
                    self.header.append(line)

        # no need for stationary point steps if not an optimisation
        if self.input.contrl.runtyp != "optimize":
            for index, line in enumerate(self.header):
                if line.split()[0] == "$STATPT":
                    del self.header[index]

        self.header = "".join(self.header)

    def change_charge_and_mult(self):
        """
        Changes charge and multiplicity unless user defines values. 
        In that case, the user-defined charge and multiplicity are used.
        """
        user_assigned_charge = False
        user_assigned_mult = False
        if hasattr(self, "user_settings"):
            user_assigned_charge = hasattr(self.user_settings, "input.contrl.icharg")
            user_assigned_mult = hasattr(self.user_settings, "input.contrl.mult")
        if self.mol.overall_charge != 0 and not user_assigned_charge:
            self.input.contrl.icharg = self.mol.overall_charge
        if self.mol.overall_mult != 1 and not user_assigned_mult:
            self.input.contrl.mult = self.mol.overall_mult

    def make_automatic_changes(self):
        """
        Common scenarios are implemented to remove the commands needed to be called by the user.
        For example, this sets the opposite spin parameter for SRS-MP2 for commonly used correlation
        consistent basis sets without the need to define it in the user code.
        """

        self.change_charge_and_mult()

        opp_spin_params = {
            "cct": 1.64,
            "ccq": 1.689,
            "accd": 1.372,
            "acct": 1.443,
            "accq": 1.591,
        }
        if "mp2" in self.input:
            for basis, opp in opp_spin_params.items():
                if self.input.basis.gbasis.lower() == basis:
                    self.input.mp2.scsopo = opp
                    break
        if self.fmo and "pcm" in self.input:
            self.input.pcm.ifmo = -1

    def parse_settings(self):
        """Transforms all contents of |Settings| objects into GAMESS input file headers, containing all the information pertinent to the calculation"""

        def format_line_if_too_long(line):
            if len(line) > 55:
                line = line.split()
                # find suitable character length to split on
                line_length = 0
                chars = []
                for val in line:
                    line_length += len(val)
                    chars.append(line_length)

                # For super long lines, could call this recursively
                for index, length in enumerate(chars):
                    if length > 55:
                        line.insert(index - 1, "\n ")
                        break  # only insert one newline, and indent
                # space before, newline at end
                line.insert(0, "")
                line = " ".join(line) + "\n"
            return line

        def preserve_value(value):
            """
            If user inserts newline, don't remove it.
            Presumably there is some desired formatting of that group
            """
            return "\n" in value

        def parse(key, value):
            ret = ""

            if isinstance(value, Settings):
                ret += " ${}".format(key.upper())
                for el in value:
                    ret += " {}={}".format(el.upper(), str(value[el]).upper())
                ret += " $END\n"
            else:
                ret += " ${} {}\n $END\n".format(key.upper(), value.upper())
            if not preserve_value(value):
                ret = format_line_if_too_long(ret)
            return ret

        inp = [parse(item, self.input[item]) for item in self.input]
        return inp

    def fmo_meta(self):
        """Creates strings for the INDAT and ICHARG blocks of GAMESS FMO calculations, bound to the
        molecule instance as self.indat and self.charg"""

        info = {}
        # group together indat and charge for each fragment, and order according
        # to the atom indices of indat
        for frag, data in self.mol.fragments.items():
            if frag is not "ionic":
                if len(data["atoms"]) == 1:
                    # should add to next fragment- wasteful to run on own node
                    info[frag] = {
                        "indat": f"0,{data['atoms'][0].index},-{data['atoms'][0].index},",
                        "charg": str(data["charge"]),
                        "mult": str(data["multiplicity"]),
                    }
                else:
                    # check for consecutive numbers
                    atom_indices = [atom.index for atom in data["atoms"]]
                    if consecutive(atom_indices):
                        indat_string = (
                            f"0,{data['atoms'][0].index},-{data['atoms'][-1].index},"
                        )
                    else:
                        # odd ordering, especially with fragmenting on a bond
                        groups = []
                        frags = []
                        indices = [atom.index for atom in data["atoms"]]
                        for i, atom_i in enumerate(indices):
                            for j, atom_j in enumerate(indices):
                                if j == i + 1:
                                    # first item
                                    if len(frags) == 0:
                                        frags.append(atom_i)
                                    if not consecutive([atom_i, atom_j]):
                                        if frags[-1] == "-":
                                            # end one fragment, start another
                                            frags += [atom_i, atom_j]
                                    # count back to see how many atoms are in 'current' fragment
                                    count = 0
                                    for val in frags[::-1]:
                                        if val == "-":
                                            count += 1
                                        else:
                                            break
                                    if count == 0:
                                        frags.append("-")
                                    # at end
                                    if j == len(indices) - 1:
                                        frags.append(atom_j)

                        # Algorithm above needs improving- if it worked properly,
                        # there would be no need for the next two for loops

                        # Check for 'single' atom parts i.e. 43, -, 43,
                        # remove -, 43
                        for i, val in enumerate(frags):
                            if frags[i : i + 3] == [val, "-", val]:
                                del frags[i + 1 : i + 3]
                        # Check for consecutives like 46, -, 47, remove - from middle
                        for i, val in enumerate(frags):
                            try:
                                if (
                                    isinstance(frags[i], int)
                                    and isinstance(frags[i + 1], str)
                                    and isinstance(frags[i + 2], int)
                                ):
                                    if frags[i + 2] == val + 1:
                                        del frags[i + 1]
                            except IndexError:
                                break
                        groups.append(frags)

                        indat_string = "0,"
                        for group in groups:
                            if len(group) == 1:
                                indat_string += f"{group[0]},"
                            else:
                                # indat_string += f'{group[0]},-{group[-1]},'
                                # if '-', then add to next element
                                for i, val in enumerate(group):
                                    if val == "-":
                                        group[i + 1] = f"-{group[i+1]}"
                                        del group[i]
                                group = map(str, group)
                                indat_string += ",".join(group) + ","

                    info[frag] = {
                        "indat": indat_string,
                        "charg": str(data["charge"]),
                        "mult": str(data["multiplicity"]),
                    }

        # items need sorting
        # 0,1,7, ### sort on 2nd item ###
        # 0,8,28,
        # 0,29,35,
        # and not
        # 0,1,7,
        # 0,29,35,
        # 0.8,28

        sorted_info = sorted(
            info.items(), key=lambda val: int(val[1]["indat"].split(",")[1])
        )
        # could also just sort on mol (or frag of info[frag]), from the assignments in self.split(), but these might not
        # always be in a numerical order- by using the index from self.coords, it is always ensured that the
        # correct order is shown, as these coords are also used in the input file
        self.fmo_indat = []
        self.fmo_charg = []
        self.fmo_mult = []
        for val in sorted_info:
            self.fmo_indat.append(val[1]["indat"])
            self.fmo_charg.append(val[1]["charg"])
            self.fmo_mult.append(val[1]["mult"])
        self.fmo_indat.append("0")

    def fmo_formatting(self):
        self.fmo_meta()  # gives self.mol.indat, self.mol.charg
        if self.input.contrl.runtyp.lower() in (
            "optimize",
            "hessian",
            "fmohess",
            "sadpoint",
        ):
            nbody = 2
            rcorsd = 100
        else:
            # FMO3 for specs- change rcorsd to 50
            nbody = 3
            rcorsd = 50

        string = f"\n     NFRAG={len(self.mol.fragments)} NBODY={nbody}\n"
        if "mp2" in self.input:
            string += "     MPLEVL(1)=2\n"
        string += f"     INDAT(1)={self.fmo_indat[0]}\n"
        for d in self.fmo_indat[1:]:
            string += f"{' '*14}{d}\n"
        string += f"     ICHARG(1)={','.join(self.fmo_charg)}\n"
        string += f"     MULT(1)={','.join(self.fmo_mult)}\n"
        string += f"     RESPAP=0 RESPPC=-1 RESDIM=100 RCORSD={rcorsd}"
        return string

    def make_inp(self):
        inp = self.header
        inp += " $DATA\n"
        inp += f"{self.title}\n"
        inp += "C1\n"
        if self.fmo:
            # list of tuples [('H', 1.0), ('O', 8.0)]
            for el in self.mol.complex["elements"]:
                inp += f" {el[0]} {el[1]}\n"
            inp += " $END\n"
            inp += " $FMOXYZ\n"
        for atom in self.mol.coords:
            inp += f" {atom.symbol:5s} {PT.get_atnum(atom):>3}.0{atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}\n"
        inp += " $END"
        return inp

    def file_basename(self):
        """If no filename is passed when the class is instantiated, the name of the file defaults to
        the run type: a geometry optimisation (opt), single point energy calculation (spec), 
        or a hessian matrix calculation for vibrational frequencies (hess). 
        This method creates an attribute ``base_name``, used in creating the input and job files."""

        if self.filename is not None:
            self.base_name = self.filename
        else:
            options = {
                "optimize": "opt",
                "energy": "spec",
                "hessian": "hess",
                "fmohess": "hess",
                "sadpoint": "ts",
            }
            self.base_name = options.get(
                self.input.contrl.runtyp, "file"
            )  # default name = file

    def create_inp(self):
        self.input = self.input.remove_none_values()
        self.determine_fragments()  # add fmo info to input settings, if self.fmo is True
        self.make_automatic_changes()
        self.unordered_header = self.parse_settings()
        self.order_header()  # create self.header variable
        inp = self.make_inp()
        self.file_basename()
        self.write_file(inp, filetype="inp")

    def get_job_template(self):
        if "dfttyp" in self.input.contrl.keys():
            job_file = self.find_job(dft=True)
        else:
            job_file = self.find_job()
        with open(job_file) as f:
            job = f.read()
            return job

    def change_mgs_job(self, job):
        if hasattr(self.mol, "fragments") and len(self.mol.fragments) != 0:
            num_frags = len(self.mol.fragments)
            jobfile = job.replace("nodes=1", f"nodes={num_frags}")
            jobfile = jobfile.replace("24 24", f"{24 * num_frags} 24")
            return jobfile
        return job

    def change_rjn_job(self, job):
        if hasattr(self.mol, "fragments") and len(self.mol.fragments) != 0:
            num_frags = len(self.mol.fragments)
            jobfile = job.replace("ncpus=32", f"ncpus={16 * num_frags}")
            jobfile = jobfile.replace(
                "mem=125gb", f"mem={4 * 16 * num_frags}gb"
            )  # 4gb cpus
            jobfile = jobfile.replace(
                "jobfs=150gb", f"jobfs={4 * 16 * num_frags + 20}gb"
            )
            return jobfile
        return job

    def change_stm_job(self, job):
        jobfile = job.replace("name", f"{self.base_name}")
        if hasattr(self.mol, "fragments") and len(self.mol.fragments) != 0:
            num_frags = len(self.mol.fragments)
            jobfile = jobfile.replace("-N 1", f"-N {num_frags}")
            jobfile = jobfile.replace("-n 22", f"-n {22 * num_frags}")
        return jobfile

    def create_job(self):
        """Returns the relevant job template as a list, then performs the necessary modifications. After, the job file is printed in the appropriate directory."""
        jobfile = self.get_job_template()
        if self.sc == "mgs":
            jobfile = self.change_mgs_job(jobfile)
            jobfile = jobfile.replace("name", f"{self.base_name}")
        elif self.sc == "rjn":
            jobfile = self.change_rjn_job(jobfile)
            jobfile = jobfile.replace("name", f"{self.base_name}")
        elif self.sc == "mon":
            jobfile = jobfile.replace("base_name", f"{self.base_name}")
        elif self.sc == "mas":
            jobfile = jobfile.replace("base_name", f"{self.base_name}")
        elif self.sc == "stm":
            jobfile = self.change_stm_job(jobfile)
        self.write_file(jobfile, filetype="job")

    def make_run_dir(self):
        if not self.made_run_dir:  # only do it once
            if not exists(self.base_name):
                mkdir(self.base_name)  # make opt/spec/hessin parent dir
            self.made_run_dir = True

    def place_files_in_dir(self):
        """Move input and job files into a directory named with the input name (``base_name``) i.e.
        moves opt.inp and opt.job into a directory called ``opt``."""
        complex_dir = join(getcwd(), "complex")
        if self.is_complex:
            if not exists(complex_dir):
                mkdir(complex_dir)
            system(f"mv {self.base_name}.inp {self.base_name}.job complex/")
            system(f"cp {self.xyz} complex/complex.xyz")

    def ionic_mol_has_two_or_more_frags(self):
        """
        If ionic molecule has more than two fragments, return True, and use fmo,
        else make a non-fmo calculation
        """
        ionic_mol = Molecule(atoms=self.mol.ionic["atoms"])
        ionic_mol.separate()
        return len(ionic_mol.fragments) > 2

    def create_inputs_for_fragments(self, complex_is_fmo=False):
        """Very useful to generate files for each fragment automatically, 
        for single point and frequency calculations, generating free energy changes. 
        Called if ``frags_in_subdir`` is set to True, as each fragment is given a 
        subdirectory in an overall subdirectory, creating the following directory 
        structure (here for a 5-molecule system):
            .
            ├── frags
            │   ├── acetate0
            │   │   ├── acetate0.xyz
            │   │   └── spec.inp
            │   ├── acetate1
            │   │   ├── acetate1.xyz
            │   │   └── spec.inp
            │   ├── choline2
            │   │   ├── choline2.xyz
            │   │   └── spec.inp
            │   ├── choline3
            │   │   ├── choline3.xyz
            │   │   └── spec.inp
            │   └── water4
            │       ├── spec.inp
            │       └── water4.xyz
            ├── spec.inp
        """
        self.is_complex = False
        # look over self.mol.fragments, generate inputs- make a settings object with the desired features
        if not hasattr(self.mol, "fragments"):
            self.mol.separate()
        # make subdir if not already there
        subdirectory = join(getcwd(), "frags")
        if not exists(subdirectory):
            mkdir(subdirectory)

        parent_dir = getcwd()
        count = 0  # avoid  overwriting files by iterating with a number
        for frag, data in self.mol.fragments.items():
            if data["frag_type"] == "frag":

                # make a directory inside the subdir for each fragment
                # i.e. acetate0, acetate1, choline2, choline3, water4
                name = f"{data['name']}_{count}"
                if not exists(join(subdirectory, name)):
                    mkdir(join(subdirectory, name))  # ./frags/water4/
                chdir(join(subdirectory, name))
                write_xyz(atoms=data["atoms"], filename=name + str(".xyz"))

                # re-use settings from complex
                if hasattr(self, "merged"):
                    frag_settings = self.merged
                else:
                    frag_settings = self.defaults
                frag_settings.input.contrl.icharg = data["charge"]
                if data["multiplicity"] != 1:
                    frag_settings.input.contrl.mult = data["multiplicity"]
                job = GamessJob(
                    using=name + str(".xyz"), settings=frag_settings, run_dir=True
                )
                chdir(parent_dir)
                count += 1

        if hasattr(self.mol, "ionic"):
            if len(self.mol.ionic["atoms"]) > 0:
                # only 1 ionic network
                subdir_ionic = join(getcwd(), "ionic")
                if not exists(subdir_ionic):
                    mkdir(subdir_ionic)
                chdir(subdir_ionic)
                write_xyz(atoms=self.mol.ionic["atoms"], filename="ionic.xyz")

                # re-use settings from complex
                if hasattr(self, "merged"):
                    frag_settings = self.merged
                else:
                    frag_settings = self.defaults
                frag_settings.input.contrl.icharg = self.mol.ionic["charge"]
                if self.mol.ionic["multiplicity"] != 1:
                    frag_settings.input.contrl.mult = self.mol.ionic["multiplicity"]

                # FMO only if more than 2 fragments
                if complex_is_fmo:
                    complex_is_fmo = self.ionic_mol_has_two_or_more_frags()

                job = GamessJob(
                    using="ionic.xyz",
                    settings=frag_settings,
                    fmo=complex_is_fmo,
                    run_dir=True,
                )
                chdir(parent_dir)
