from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.settings import Settings, read_template, dict_to_settings
from ..core.job import Job
from ..core.periodic_table import PeriodicTable as PT
from ..core.sc import Supercomp
from ..core.utils import search_dict_recursively, write_xyz

from os import chdir, mkdir, getcwd, system
from os.path import exists, join, dirname

__all__ = ["PsiJob"]


class PsiJob(Job):
    # Note the job scripts require the supercomputer to be entered, such as:

    # >>> j = PsiJob(using = 'file.xyz')
    # >>> j.supercomp = 'raijin'
    """Class for creating PSI4 input files and job scripts. 
    
    The input files generated default to single point energy calculations using MP2/cc-pVTZ, with frozen core orbitals- for this reason, these single point calculations are very fast. This is easily changed by creating a |Settings| object and adding parameters, using the following syntax.
    >>> s = Settings()
    >>> s.input.globals.basis= 'cc-pVDZ'
    >>> s.input.molecule.extra_value = 'extra'
    >>> j = PsiJob(using = '../xyz_files/mesylate.xyz', settings = s)
    This yields the following result:
        memory 32 Gb
        molecule complex {
        -1 1
         C      -11.52615475      2.13587901     -3.92614475
         H      -12.17727298      1.37283268     -4.39314733
         H      -12.13111156      2.84527650     -3.33020803
         H      -10.95289836      2.67720525     -4.70258006
         S      -10.36648767      1.31567304     -2.82897636
         O       -9.54405868      2.38757303     -2.22205822
         O      -11.24567273      0.60890457     -1.83183396
         O       -9.60100212      0.36690604     -3.68579623
        units angstrom
        no_reorient
        symmetry c1
        extra_value extra
        }
        set globals {
            basis cc-pVDZ
            scf_type DF
            freeze_core True
            guess sad
            S_ORTHOGONALIZATION canonical
        }
        energy('mp2')
    Options are added in sections: 
        - self.input.molecule for any key value pair in the molecule section
        - self.input.unbound for any key value pair outside of molecule and globals. The value can be a string or a list.
            >>> self.input.unbound.key = 'value'
            # key value
             >>> self.input.unbound.key2 = 'value value value'
            # key value value value
            >>> self.input.unbound.key = ['value1', 'value2', 'value3']
            # key value1 value2 value3
        - self.input.globals for the 'set globals' section
        - any options not enclosed in braces appear before the last line
        - To change the run type:
            >>> self.input.run = {'optimize': 'scf'}
            # optimize('scf')
        - If extra run options are required:
            >>> self.input.run.additional = {'dertype': 'energy'} 
            # optimize('scf', dertype='energy')
            >>> self.input.run.additional = {'dertype': 'energy', 'option2': 'value'} 
            # optimize('scf', dertype='energy', 'option2'='value')
        NOTE: An option for adding commands outside of the molecule and globals section needs to be added.                        
    
    The names of files created default to the type of calculation: optimisation (opt), single point
energy (spec) or hessian matrix calculation for thermochemical data and vibrational frequencies (hess). If a different name is desired, pass a string with the ``filename`` parameter, with no extension. The name will be used for both input and job files.
        >>> job = PsiJob(using = 'file.xyz', filename = 'benzene')
    This command produces two files, benzene.inp and benzene.job.
    
    To run a counterpoise corrected calculation, pass in `cp=True` to the constructor:
        >>> job = PsiJob(using = 'file.xyz', filename = 'benzene', cp=True)
    
   This produces an extra file- a counterpoise corrected Hartree-Fock calculation of the entire
cluster.
       
 
    """

    def __init__(
        self,
        using=None,
        frags_in_subdir=False,
        settings=None,
        filename=None,
        is_complex=False,
        run_dir=None,
        cp=False,
    ):
        super().__init__(using)
        self.filename = filename
        self.defaults = read_template("psi.json")  # settings object
        if settings is not None:
            # can only have one run type, currently- need to delete the energy if running an
            # optimisation, for example
            if "run" in settings.input.keys():
                del self.defaults.input.run
            self.merged = self.defaults.merge(settings)  # merges inp, job data
            self.input = self.merged.input
            self.job = self.merged.job
            self.meta = self.merged.meta
        else:
            self.input = self.defaults.input
            self.meta = self.defaults.meta  # just to save hasattr(self, 'meta') further down
        if "/" in using:
            self.title = using.split("/")[-1][:-4]  # say using = ../xyz_files/file.xyz -->
        else:
            self.title = using[:-4]

        self.is_complex = is_complex  # creates a `complex` dir

        self.cp = cp
        self.create_inp(counterpoise=self.cp)
        self.create_job()
        self.place_files_in_dir()
        if frags_in_subdir:
            self.create_inputs_for_fragments()

    def make_header(self):
        """Transform all contents of |Settings| objects into PSI4 input file headers, containing all the information pertinent to the calculation"""
        self.find_charge_and_mult()
        comment = f"# PSI4 Calc: {self.title}\n\n"
        mem = f"memory {self.input.memory}\n\n"
        mol = "molecule complex {\n"
        charge = f"{self.input.charge} {self.input.mult}\n"
        atoms = ""
        for atom in self.mol.coords:
            atoms += f" {atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}\n"
        units = f"units {self.input.molecule.units}\n"
        sym = f"symmetry {self.input.molecule.symmetry}\n"
        reorient = "no_reorient\n"
        end = "}\n"

        data = [comment, mem, mol, charge, atoms, units, reorient, sym, end]

        # rm unneccesary options
        self.input = self.input.remove_none_values()

        # add in user options
        for key, value in self.input.molecule.items():
            if key not in ("charge", "multiplicity", "units", "symmetry"):
                key = f"{key} {value}\n"
                data.insert(-1, key)  # insert before last item
        self.inp = data

    def add_unbound(self):
        """May never be required- but this adds options between the molecule and global sections.
        Returns a dictionary of terms- might need more than two terms on same line = nested dict """

        vals = search_dict_recursively(self.input.unbound)
        if vals != {}:  # if not empty
            self.inp.append("\n")
            for key, value in vals.items():
                if isinstance(value, list):
                    self.inp.append(f"{key} {' '.join(value)}\n")
                elif isinstance(value, str):
                    self.inp.append(f"{key} {value}\n")

    def add_globals(self):
        self.inp.append("\nset globals {\n")
        for key, value in self.input.globals.items():
            self.inp.append(f"    {key} {value}\n")
        self.inp.append("}\n")

    def add_run(self):
        res = []
        # list of tuples- to ensure the 'normal' entry, the one defined input.run, appears first
        # in the list by adding a counter. In testing, {optimize: scf} with additional {dertype:
        # energy} produced dertype('energy', optimize='scf'), not optimize('scf', dertype='energy')
        # due to alphabetical ordering of [('dertype', 'energy'), ('optimize', 'scf')]
        # probably should have just made two lists of tuples, one for normal, one for additional
        for k, v in self.input.run.items():
            if k != "additional":
                res.append((0, k, v))
            if k == "additional":
                counter = 1
                for k1, v1 in self.input.run[k].items():
                    res.append((counter, k1, v1))
                    counter += 1
                # if I ever need to add two different types of run in the same file,
                # comment out the two lines above, add in the 3 below, and add an additional
                # dict key: (AND CHANGE DOCSTRING)
                # would need to add resulting string to a list and then concatenate that list with
                # self.inp as well, otherwise you would have a combination in the same line
                ##############
                # s = Settings()
                # s.input.run = {'optimize': 'scf'}
                # s.input.run.additional = {'optimize' :{'dertype': 'energy',
                #                                        'entry': 'value'}}
                ##############
                # for data in self.input.run[k].values():
                #     for k1, v1 in data.items():
                #         res.append((k1, v1))
        res = sorted(res, key=lambda val: val[0])  # sort by the first item of tuple, the number
        string = f"{res[0][1]}('{res[0][2]}'"
        for val in res[1:]:
            string += f", {val[1]}='{val[2]}'"
        string += ")"
        self.inp.append(string)
        self.inp = "".join(self.inp)

    def file_basename(self):
        """If no filename is passed when the class is instantiated, the name of the file defaults to
        the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (hess). This method creates an attribute ``base_name``, used in creating the input and job files."""
        for key in self.input.run.keys():  # run, or additional
            if key != "additional":
                nom = key
        if self.filename == None:
            options = {"optimize": "opt", "energy": "spec", "frequency": "hess"}
            self.base_name = options.get(nom, "file")  # default name = file
        else:
            self.base_name = self.filename

    def write_file(self, data, filetype):
        """Writes the generated PSI4 input/jobs to a file. If no filename is passed when the class is instantiated, the name of the file defaults to the run type: a geometry optimisation (opt), single point energy calculation (spec), or a hessian matrix calculation for vibrational frequencies (freq). 
        NOTE: Must pass data as a string, not a list!"""
        with open(f"{self.base_name}.{filetype}", "w") as f:
            f.write(data)

    def make_counterpoise(self):
        """
        Make a counterpoise corrected HF input file and place in a separate directory.
        """
        # make new PsiJob object, the only thing that changes is atoms, and directory name
        # class PsiJob_CP(PsiJob):
        #
        #   def convert_atoms_to_cp(self):
        #
        #       if not hasattr(self.mol, 'fragments'):
        #           self.mol.separate()
        #       for frag in self.mol.fragments.values():
        #           for atom in frag['atoms']:
        #               atoms.append(f" {atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}")
        #           atoms.append('--')
        #       atoms = "\n".join(atoms[:-1]) + '\n'
        #       return atoms

        self.find_charge_and_mult()
        comment = f"# PSI4 Calc: {self.title}\n\n"
        mem = f"memory {self.input.memory}\n\n"
        mol = "molecule complex {\n"
        charge = f"{self.input.charge} {self.input.mult}\n"
        atoms = []
        if not hasattr(self.mol, "fragments"):
            self.mol.separate()
        for frag in self.mol.fragments.values():
            for atom in frag["atoms"]:
                atoms.append(f" {atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}")
            atoms.append("--")
        atoms = "\n".join(atoms[:-1]) + "\n"
        units = f"units {self.input.molecule.units}\n"
        sym = f"symmetry {self.input.molecule.symmetry}\n"
        reorient = "no_reorient\n"
        end = "}\n"

        data = [comment, mem, mol, charge, atoms, units, reorient, sym, end]

        # add in user options
        for key, value in self.input.molecule.items():
            if key not in ("charge", "multiplicity", "units", "symmetry"):
                key = f"{key} {value}\n"
                data.insert(-1, key)  # insert before last item
        data.append("\nset globals {\n")
        for key, value in self.input.globals.items():
            data.append(f"    {key} {value}\n")
        data.append("}\n")

        data.append("energy('HF', bsse_type='cp')")
        cp_dir = join(getcwd(), "cp-hf")
        if not exists(cp_dir):
            mkdir(cp_dir)
        cp_input = join(cp_dir, f"{self.base_name}.inp")
        cp_job = join(cp_dir, f"{self.base_name}.job")
        with open(cp_input, "w") as f:
            for line in data:
                f.write(line)

        job_file = self.find_job()
        with open(job_file) as f:
            job = f.read()

        if self.sc == "mgs":
            job = job.replace("name", f"{self.base_name}")
        elif self.sc == "rjn":
            # should alter the job time as they never need 4 hours-
            # walltime = max_time_for_4ip (probs have?) * num atoms / num atoms in 4IP
            job = job.replace("name", f"{self.base_name}")
        elif self.sc == "mas":
            job = job.replace("base_name", f"{self.base_name}")
        elif self.sc == "mon":
            job = job.replace("base_name", f"{self.base_name}")
        elif self.sc == "stm":
            job = job.replace("name", f"{self.base_name}")

        if "time" in self.meta:
            job = job.replace("24:00:00", self.meta.time)

        # pbs, also needs adding in normal job file and for slurm
        if self.sc in super().PBS_HOSTS:
            if "nproc" in self.meta:
                job = job.replace("ncpus=16", f"ncpus={self.meta.nproc}")
            if "jobfs" in self.meta:
                jobfs = self.meta.jobfs[:-2]  # drop units
                job = job.replace("jobfs=10GB", f"jobfs={jobfs}GB")
            if "mem" in self.meta:
                mem = self.meta.mem[:-2]
                job = job.replace("mem=64GB", f"mem={mem}GB")

        with open(cp_job, "w") as j:
            j.write(job)

    def create_inp(self, counterpoise=False):
        self.make_header()
        self.add_unbound()
        self.add_globals()
        self.add_run()
        self.file_basename()
        self.write_file(self.inp, filetype="inp")
        if counterpoise:
            self.make_counterpoise()

    def get_job_template(self):
        job_file = self.find_job()
        with open(job_file) as f:
            job = f.read()
            return job

    def create_job(self):
        """Returns the relevant job template as a list, then performs the necessary modifications. After, the job file is printed in the appropriate directory."""
        jobfile = self.get_job_template()
        # modify
        if self.sc == "mgs":
            jobfile = jobfile.replace("name", f"{self.base_name}")
        elif self.sc == "rjn":
            # should alter the job time as they never need 4 hours-
            # walltime = max_time_for_4ip (probs have?) * num atoms / num atoms in 4IP
            jobfile = jobfile.replace("name", f"{self.base_name}")
        elif self.sc == "mas":
            jobfile = jobfile.replace("base_name", f"{self.base_name}")
        elif self.sc == "mon":
            jobfile = jobfile.replace("base_name", f"{self.base_name}")
        elif self.sc == "stm":
            jobfile = jobfile.replace("name", f"{self.base_name}")

        if "time" in self.meta:
            job = job.replace("24:00:00", self.meta.time)

        self.write_file(jobfile, filetype="job")

    def place_files_in_dir(self):
        """
        Move input and job files into a directory named `complex`, if self.is_complex
        is set to True
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
        count = 0  # avoid  overwriting files by iterating with a number
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
                frag_settings.input.molecule.charge = data["charge"]
                if data["multiplicity"] != 1:
                    frag_settings.input.molecule.multiplicity = data["multiplicity"]
                job = PsiJob(using=name + str(".xyz"), settings=frag_settings, run_dir=True)
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
            frag_settings.input.molecule.charge = self.mol.ionic["charge"]
            if self.mol.ionic["multiplicity"] != 1:
                frag_settings.input.molecule.multiplicity = self.mol.ionic["multiplicity"]
            job = PsiJob(using="ionic.xyz", settings=frag_settings, run_dir=True)
            chdir(parent_dir)
