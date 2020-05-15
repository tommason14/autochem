from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.settings import Settings, read_template
from ..core.job import Job
from os import mkdir, chdir, getcwd
from os.path import exists, join
from shutil import copyfile, move

__all__ = ["GaussJob"]


class GaussJob(Job):
    """Class for creating Gaussian input files and job scripts. 
    
    The names of files created default to the type of calculation: 
    optimisation (opt), single point energy (spec) or hessian matrix 
    calculation for thermochemical data and vibrational frequencies (hess). 
    If a different name is desired, pass a string with the ``filename`` 
    parameter, with no extension. The name will be used for both input and job
    files.
    
        >>> job = GaussJob(using = 'file.xyz', filename = 'benzene')
    
    This command produces 'benzene.job', containing both input data and 
    job scheduler information.
    For meta data, number of processors (ncpus), memory (mem) etc, use sett.meta.ncpus=46.
    """

    def __init__(self, using=None, frags_in_subdir=False, settings=None, filename=None):
        super().__init__(using)
        self.filename = filename
        self.defaults = read_template("gaussian.json")
        if settings is not None:
            self.user_settings = settings.as_dict()
            self.merged = self.defaults.merge(settings)
            self.meta = self.merged.meta
            self.input = self.merged.input
            self.frag = Settings()
            self.frag.meta = self.merged.frag.meta
        else:
            self.input = self.defaults.input
            self.meta = self.defaults.meta
            self.frag = Settings()
            self.frag.meta = self.merged.frag.meta
        self.input = self.input.remove_none_values()
        if "/" in using:
            self.title = using.split("/")[-1][:-4]
        else:
            self.title = using[:-4]
        self.xyz = using

        self.file_basename()
        if frags_in_subdir:
            self.create_inputs_for_fragments()
        self.write_file(self.inp, filetype="job")

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
            if self.runtype == "":
                self.base_name = "spec"
            elif "opt" in self.runtype and "freq" in self.runtype:
                self.base_name = "opt-freq"
            elif "opt" in self.runtype and "freq" not in self.runtype:
                self.base_name = "opt"
            else:
                self.base_name = "freq"

    @property
    def job_data(self):
        def _change_partition(jobfile, partition, search_term):
            for num, line in enumerate(jobfile):
                if search_term in line:
                    jobfile[num] = f"{search_term}{partition}"
            return jobfile

        job = self.get_job_template().replace("name", self.base_name)
        if hasattr(self, "meta"):
            if "time" in self.meta:
                job = job.replace("24:00:00", self.meta.time)
            if "nodemem" in self.meta:
                mem = self.meta.nodemem[:-2]
                if self.sc in ("mas", "mon"):
                    job = job.replace("mem=32", f"mem={mem}")
                if self.sc == "gadi":
                    job = job.replace("mem=192", f"mem={mem}")
            if "ncpus" in self.meta:
                if self.sc in super().SLURM_HOSTS:
                    job = job.replace(
                        "cpus-per-task=16", f"cpus-per-task={self.meta.ncpus}"
                    )
                    # for stampede, specified as -c, so it won't change there, which is
                    # what we want as you are charged for the whole node there!
                else:  # gadi
                    job = job.replace("ncpus=48", f"ncpus={self.meta.ncpus}")
            if "partition" in self.meta:
                if self.sc in super().SLURM_HOSTS:
                    jobfile = job.split("\n")
                    jobfile = _change_partition(
                        jobfile, self.meta.partition, search_term="#SBATCH -p "
                    )
                    # might be --partition=
                    jobfile = _change_partition(
                        jobfile, self.meta.partition, search_term="#SBATCH --partition="
                    )
                    job = "\n".join(jobfile)
                else:
                    job = job.replace(
                        "#PBS -l wd", f"#PBS -q {self.meta.partition}\n#PBS -l wd"
                    )

            if self.sc in super().PBS_HOSTS:
                if "jobfs" in self.meta:
                    jobfs = self.meta.jobfs.upper().replace("GB", "")
                    job = job.replace("jobfs=200GB", f"jobfs={jobfs}GB")

        return job

    @property
    def inp(self):
        """
        Creates Gaussian input file of the form:
        SLURM/PBS data
        <blank>
        name
        <blank>
        charge/mult
        xyzdata
        <blank>
        END
        """
        inp = [
            self.job_data,
            self.metadata,
            self.run_info,
            self.title,
            self.coord_info,
            "END",
        ]
        return "\n\n".join(inp)

    @property
    def metadata(self):
        """
        Include data such as memory and number of cpus in the Gaussian file.
        """
        excluded_properties = ("time", "partition", "nodemem", "jobfs")
        # input by user for scheduler
        meta = []
        if self.sc == "stm":
            self.meta.mem = "160gb"
            self.meta.ncpus = 46
        for k, v in self.meta.items():
            if k not in excluded_properties:
                if k == "ncpus":
                    meta.append(f"%nproc={v}")
                else:
                    meta.append(f"%{k}={v}")

        return "\n".join(meta).replace("name", self.base_name)

    @property
    def runtype(self):
        """
        Decides if the job should be an optimisation, single point,
        frequency job, or an optimisation followed by a frequency job.
        Also adds parameters like ts,eigentest to the opt call, for example.
        """
        params = self.input.keys()
        runtype = ""
        if "opt" in params:
            runtype += gauss_print(self.input, "opt")
        if "freq" in params:
            if len(runtype) > 0:
                runtype += " "
            runtype += "freq"
        if "opt" not in params and "freq" not in params:
            runtype = ""
        return runtype

    @property
    def additional_params(self):
        """
        Add in parameters to the `run_info` that do not involve 
        a basis set, method or run type.
        """
        addn = ""
        ignore = ("opt", "freq", "method", "basis", "meta", "charge", "mult")
        for arg in self.input:
            if arg not in ignore:
                addn += gauss_print(self.input, arg) + " "
        return addn

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
    def run_info(self):
        """
        Creating a line like this:
        #P wB97XD/cc-pVDZ opt=(ts,noeigentest,calcfc) freq SCF=tight SCRF=(SMD,solvent=water) INT=(grid=ultrafine)
        from data stored in self.input
        """
        return (
            f"#P {self.input.method}/{self.input.basis}"
            f"{self.formatted_run}{self.additional_params}"
        )

    @property
    def coord_info(self):
        self.find_charge_and_mult()
        info = [f"{self.input.charge} {self.input.mult}"]
        info += [
            f"{atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}"
            for atom in self.mol.coords
        ]
        return "\n".join(info)

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
                # for job info, use self.frag.meta
                frag_settings = frag_settings.merge(self.frag)
                frag_settings.input.charge = data["charge"]
                if data["multiplicity"] != 1:
                    frag_settings.input.mult = data["multiplicity"]
                job = GaussJob(using=name + str(".xyz"), settings=frag_settings)
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
            job = GaussJob(using="ionic.xyz", settings=frag_settings)
            chdir(parent_dir)


def gauss_print(d, value):
    """
    Decides how to print a parameter.
    For example, opt, opt=ts or opt=(ts,eigentest,calcfc). 
    Checks a |Settings| object, `d`, for a key, `value`.
    For example, sett.input.freq=True in the script would
    produce 'freq', sett.input.scf='tight' would produce
    'scf=tight', and sett.input.opt='ts,noeigentest,calcfc' 
    produces 'opt=(ts,noeigentest,calcfc)'. Note: doesn't
    work with lists or dict values, but unlikely that they would
    be passed in as settings values anyway.
    """
    if isinstance(d[f"{value}"], bool) or isinstance(d[f"{value}"], int):
        return value
    elif len(d[f"{value}"].split(",")) > 1 or "=" in d[f"{value}"]:
        return f"{value}=({d[f'{value}']})"
    else:
        return f"{value}={d[f'{value}']}"
