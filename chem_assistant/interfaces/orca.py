from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.settings import Settings, read_template, dict_to_settings
from ..core.job import Job
from ..core.periodic_table import PeriodicTable as PT
from ..core.sc import Supercomp
from ..core.utils import search_dict_recursively, write_xyz

from os import chdir, mkdir, getcwd, system
from os.path import exists, join, dirname
import re

__all__ = ["OrcaJob"]


class OrcaJob(Job):
    """Class for creating Orca input files and job scripts. 
    
    The input files generated default to single point energy calculations, with density fitting. This is easily changed by creating a |Settings| object and adding parameters, using the following syntax.
    >>> s = Settings()
    >>> s.input.basis='cc-pVDZ'
    >>> s.input.method=
    >>> j = OrcaJob(using='mesylate.xyz', settings = s)
    For solvent effects:
    >>> s.input.solvent.model='cpcm'
    >>> s.input.solvent.dielectric=80.0
    Or:
    >>> s.input.solvent.model='cpcm'
    >>> s.input.solvent.molecule='water'
    """

    _procs = {"stm": 46, "mon": 16, "mas": 16, "gadi": 46, "mgs": 24}

    def __init__(
        self,
        using=None,
        frags_in_subdir=False,
        settings=None,
        filename=None,
        is_complex=None,
    ):
        super().__init__(using)
        self.filename = filename
        self.defaults = read_template("orca.json")  # settings object
        if settings is not None:
            self.user_settings = settings.as_dict()
            self.merged = self.defaults.merge(settings)  # merges inp, job data
            self.input = self.merged.input
            self.meta = self.merged.meta
        else:
            self.input = self.defaults.input
        if "/" in using:
            self.title = using.split("/")[-1][:-4]
        else:
            self.title = using[:-4]

        self.xyzfile = using.split("/")[-1]

        self.file_basename()
        self.get_sc()  # required to be called here as func uses sett.supercomp if provided

        self.is_complex = is_complex  # creates a `complex` dir

        self.write_file(self.inp, filetype="inp")
        self.create_job()
        self.place_files_in_dir()
        if frags_in_subdir:
            self.create_inputs_for_fragments()

    def write_file(self, data, filetype):
        """Writes the generated Orca input/jobs to a file. If no filename is passed when the class is instantiated, the name of the file defaults to the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (freq). 
        NOTE: Must pass data as a string, not a list!"""
        with open(f"{self.base_name}.{filetype}", "w") as f:
            f.write(data)

    def get_job_template(self):
        job_file = self.find_job()
        with open(job_file) as f:
            job = f.read()
            return job

    def create_job(self):
        """
        Returns the relevant job template as a list, then performs the 
        necessary modifications. After, the job file is printed in the      
        appropriate directory.
        """
        jobfile = self.get_job_template().replace(" name", f" {self.base_name}")
        # need space due to `dirname` being used in script

        # replace cpus for parallelisation
        # and set cpus per task to 1
        # create list for ease of manipulation
        jobfile = jobfile.split("\n")
        for num, line in enumerate(jobfile):
            # SLURM
            if re.search("SBATCH -n(tasks)? [0-9]+", line):
                jobfile[num] = f"#SBATCH -n {OrcaJob._procs[self.sc]}"
            if re.search("SBATCH -c(pus-per-task)? [0-9]+", line):
                jobfile[num] = "#SBATCH -c 1"
            # PBS?
        jobfile = "\n".join(jobfile)

        # change job time
        if "time" in self.meta:
            jobfile = jobfile.replace("24:00:00", self.meta.time)

        self.write_file(jobfile, filetype="job")

    def place_files_in_dir(self):
        """
        Move input and job files into a directory named `complex`, if self.
        is_complex is set to True
        """
        complex_dir = join(getcwd(), "complex")
        if self.is_complex:
            if not exists(complex_dir):
                mkdir(complex_dir)
            system("cp *.xyz complex/complex.xyz")
            system(f"mv {self.base_name}.inp {self.base_name}.job complex/")

    def create_inputs_for_fragments(self):
        """Very useful to generate files for each fragment automatically, for single point and frequency calculations, generating free energy changes. Called if ``frags_in_subdir`` is set to True, as each fragment is given a subdirectory in an overall subdirectory, creating the following directory structure (here for a 5-molecule system):
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
        # not necessarily any splitting prior to this
        self.is_complex = False

        self.mol.separate()  # creating frags
        # look over self.mol.fragments, generate inputs- make a settings object with the desired features

        # after separation- create another frag with the ionic cluster!

        # make subdir if not already there
        subdirectory = join(getcwd(), "frags")
        if not exists(subdirectory):
            mkdir(subdirectory)

        parent_dir = getcwd()
        count = 0  # avoid overwriting files by iterating with a number
        for frag, data in self.mol.fragments.items():
            if data["frag_type"] == "frag":
                # make a directory inside the subdir for each fragment
                name = f"{data['name']}_{count}"  # i.e. acetate0, acetate1, choline2, choline3, water4
                if not exists(join(subdirectory, name)):
                    mkdir(join(subdirectory, name))  # ./frags/water4/
                chdir(join(subdirectory, name))
                Molecule.write_xyz(
                    self, atoms=data["atoms"], filename=name + str(".xyz")
                )  # using the method, but with no class

                # use the same settings, so if runtype is freq, generate freq inputs for all fragments too.
                if hasattr(self, "merged"):
                    frag_settings = self.merged
                else:
                    frag_settings = self.defaults
                frag_settings.input.charge = data["charge"]
                if data["multiplicity"] != 1:
                    frag_settings.input.mult = data["multiplicity"]
                job = OrcaJob(using=name + str(".xyz"), settings=frag_settings)
                chdir(parent_dir)
                count += 1
        if hasattr(self.mol, "ionic"):
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
            frag_settings.input.charge = self.mol.ionic["charge"]
            if self.mol.ionic["multiplicity"] != 1:
                frag_settings.input.mult = self.mol.ionic["multiplicity"]
            job = OrcaJob(using="ionic.xyz", settings=frag_settings)
            chdir(parent_dir)

    def file_basename(self):
        """
        If no filename is passed when the class is instantiated, the name of the
        file defaults to the run type: a geometry optimisation (opt), single
        point energy calculation (spec), or a hessian matrix calculation for
        vibrational frequencies (freq). This method creates an attribute
        ``base_name``, used in creating the input file."""

        if self.filename is not None:
            self.base_name = self.filename
        else:
            if "run" in self.input:
                if "opt" in self.input.run.lower():
                    self.base_name = "opt"
                elif "freq" in self.input.run.lower():
                    self.base_name = "freq"
                else:
                    self.base_name = "spec"
            else:
                if self.runtype == "":
                    self.base_name = "spec"
                elif "opt" in self.runtype:
                    self.base_name = "opt"
                else:
                    self.base_name = "freq"

    # def move_to_subdir(self):
    #     """
    #     Makes a subdirectory based on either the filename passed in,
    #     or the xyz file used to create the input file.
    #     The xyz file is then copied over and the input file is copied
    #     to the new directory.
    #     As a result, this method must be called after making the input file.
    #     """
    #     if self.filename is not None:
    #         newdir = join(getcwd(), self.filename)
    #     else:
    #         newdir = join(getcwd(), self.title)
    #     if not exists(newdir):
    #         mkdir(newdir)
    #     copyfile(self.xyz, newdir)
    #     move(f'{self.base_name}.job', newdir)

    @property
    def job_data(self):
        return self.get_job_template().replace("name", self.base_name)

    @property
    def inp(self):
        """
        Creates Orca input file of the form:
        ! runtype method basis density_fitting

        *xyzfile charge multiplicity xyzfile

        Adds a cpcm section if necessary.
        """
        self.input.remove_none_values()
        self.solvation_effects()
        inp = [self.run_info, self.additional_info, self.coord_info]
        if hasattr(self, "cpcm_opts"):
            inp.insert(2, self.cpcm_opts)  # after additional info
        return "\n\n".join(inp)

    @property
    def runtype(self):
        """
        Decides if the job should be an optimisation, single point or 
        frequency job.
        """
        params = self.input.keys()
        runs = ("Opt", "Freq", "NumFreq")
        for run in runs:
            if run in params:
                return run
        else:
            return ""  # spec

    def get_solvent(self):
        """
        Sees if sett.input.solvent is defined
        """
        if "solvent" in self.input.keys() and "molecule" in self.input.solvent.keys():
            return "solvent"
        else:
            return "dielectric"

    def solvation_effects(self):
        """
        Checks for sett.input.solvent, for cpcm.
        If a solvent is defined, adds `CPCM(solvent)`
        """
        solv = self.get_solvent()
        if "solvent" in self.input.keys():
            if solv == "solvent" and self.input.solvent.model.lower() == "cpcm":
                self.solvation = f"CPCM({self.input.solvent.molecule})"
            elif self.input.solvent.model.lower() == "cpcm":
                self.solvation = "CPCM"
            # add cpcm section
            extra_params = self.input.solvent.keys()
            ignore = ("model", "molecule")
            for param in ignore:
                if param in extra_params:
                    del self.input.solvent[param]
                # extra_params.remove(param)
            if len(extra_params) > 0:
                self.cpcm_opts = "%cpcm\n"
                for param in extra_params:
                    self.cpcm_opts += f"{param} {self.input.solvent[param]}\n"
                self.cpcm_opts += "end"
        else:
            self.solvation = ""

    @property
    def formatted_run(self):
        """
        If a single point calculation is desired, a runtype of '' leaves a
        double space, if '{...} {self.runtype} {...}' is used in `run_info`.
        Instead, use this function and include with no spaces in `run_info`.
        i.e. '{...}{self.formatted_run}{...}'
        """
        return " " if self.runtype == "" else f" {self.runtype} "

    @property
    def additional_params(self):
        """
        Add in parameters to the `run_info` that do not involve 
        a basis set, method or run type.
        """
        addn = ""
        ignore = (
            "runtype",
            "Opt",
            "Freq",
            "NumFreq",
            "method",
            "basis",
            "density_fitting",
            "charge",
            "mult",
            "solvent",
            "meta",
        )
        for arg in self.input:
            if arg not in ignore:
                addn += arg + " "
        return addn

    @property
    def run_info(self):
        """
        Creating a line like this:
        ! opt HF aug-cc-pVTZ RIJK additional_params
        from data stored in self.input.
        For ease, can define a `sett.input.run` string which will
        be used.
        """
        if "run" in self.input:
            return f"!{self.input.run}"
        if self.input.density_fitting is not None:
            return (
                f"!{self.formatted_run}{self.input.method} "
                f"{self.input.basis} {self.input.density_fitting} "
                f"{self.solvation} {self.additional_params}"
            )
        else:
            return (
                f"!{self.formatted_run}{self.input.method} "
                f"{self.input.basis} "
                f"{self.solvation} {self.additional_params}"
            )

    @property
    def coord_info(self):
        self.find_charge_and_mult()
        return f"*xyzfile {self.input.charge} {self.input.mult} {self.xyzfile}\n\n"

    @property
    def additional_info(self):
        """
        Add information such as:
        %pal 
            nprocs 46 
        end
        with this function.
        By default, Orca is set to run in parallel with the following number of CPUs
        per node depending on the supercomputer:
        
        - Stampede: 46        
        - Magnus: 24
        - Monarch: 16
        - Massive: 16
        - Gadi: 46

        Can also set manually with `sett.input.meta.pal='nprocs 20'`, for example.
        """
        ret = ""
        if "pal" not in self.input.meta:
            ret += f"%pal\n  nprocs {OrcaJob._procs[self.sc]}\nend"
        for key, val in self.input.meta.items():
            ret += f"\n\n%{key}\n{val}\nend"
        return ret
