from ..interfaces.gamess import GamessJob
from ..interfaces.gaussian import GaussJob
from ..interfaces.orca import OrcaJob
from ..interfaces.psi import PsiJob

import os
import glob
from shutil import copyfile

__all__ = ["xyz_to_tree"]


def get_xyz():
    return [file for file in os.listdir(".") if file.endswith(".xyz")]


def ask_package():
    print("Which software would you like?")
    print(
        """\
1. GAMESS
2. PSI4
3. GAUSSIAN
4. ORCA"""
    )
    done = False
    while not done:
        choice = int(input("Choice [1,2,3,4]: "))
        if choice in (1, 2, 3, 4):
            done = True
        else:
            print("Please choose 1-4")
    fmo = False
    if choice == 1:
        done_with_fmo = False
        while not done_with_fmo:
            run_fmo = input("Run FMO calculations? [y/n] ")
            if run_fmo.lower() in ("y", "n"):
                done_with_fmo = True
            else:
                print('Please choose "y" or "n"')
        if run_fmo == "y":
            choice = 5
    # TODO: REFACTOR CODE BELOW- DEFINITELY A BETTER WAY TO DO THIS THAN TO MAKE NEW OPTIONS, ALL IT IS
    # ACHIEVING IS SETTING FRAGS_IN_SUBDIR TO FALSE
    done_with_frags = False
    while not done_with_frags:
        frags = input("Place inputs of fragments in subdirectories? [y/n] ")
        if frags.lower() in ("y", "n"):
            done_with_frags = True
        else:
            print('Please choose "y" or "n"')
    if frags == "n":
        if choice == 1:
            choice = 6
        elif choice == 2:
            choice = 7
        elif choice == 3:
            choice = 8
        elif choice == 4:
            choice = 9
        elif choice == 5:
            choice = 10
    options = {
        1: "gamess",
        2: "psi4",
        3: "gauss",
        4: "orca",
        5: "gamess_fmo",
        6: "gamess_no_frags",
        7: "psi4_no_frags",
        8: "gauss_no_frags",
        9: "orca_no_frags",
        10: "gamess_fmo_no_frags",
    }
    return options[choice]


def job_type(package, xyz, s):

    if package == "gamess":
        return GamessJob(using=xyz, frags_in_subdir=True, settings=s, is_complex=True)
    elif package == "gamess_fmo":
        return GamessJob(
            using=xyz, fmo=True, frags_in_subdir=True, settings=s, is_complex=True
        )
    elif package == "psi4":
        return PsiJob(using=xyz, frags_in_subdir=True, settings=s, is_complex=True)
    elif package == "gauss":
        return GaussJob(using=xyz, settings=s, frags_in_subdir=True, is_complex=True)
    elif package == "orca":
        return OrcaJob(using=xyz, settings=s, frags_in_subdir=True)
    elif package == "gamess_no_frags":
        return GamessJob(using=xyz, frags_in_subdir=False, settings=s)
    elif package == "psi4_no_frags":
        return PsiJob(using=xyz, frags_in_subdir=False, settings=s)
    elif package == "gauss_no_frags":
        return GaussJob(using=xyz, settings=s, frags_in_subdir=False)
    elif package == "gamess_fmo_no_frags":
        return GamessJob(using=xyz, fmo=True, frags_in_subdir=False, settings=s)
    elif package == "orca_no_frags":
        return OrcaJob(using=xyz, settings=s, frags_in_subdir=False)


def make_dir_list(file):
    filename = file[:-4]  # rm .xyz
    new_dirs = []
    for part in filename.split("_"):
        new_dirs.append(part)
    return new_dirs


def change_to_subdir(subdirectory):
    new_dir = os.path.join(os.getcwd(), subdirectory)
    if os.path.isdir(new_dir):
        os.chdir(new_dir)
    elif not os.path.isdir(new_dir):
        os.mkdir(new_dir)
        os.chdir(new_dir)


def copy_xyz(xyz_dir, file):
    xyz = os.path.join(xyz_dir, file)
    dest = os.path.join(os.getcwd(), file)
    copyfile(xyz, dest)


def make_tree_and_copy(xyz_dir, files):
    """Makes a sibling directory to 'files' (named 'calcs'), creates subdirectories with names based on the xyz files in 'files', and then copies the xyz files from 'files' to the deepest sub directory of the path created.

    The end result is this:
    .
    ├── calcs
    │   ├── c1mim
    │   │   └── nh3
    │   │       ├── c1mim_nh3.xyz
    │   │       └── frags
    │   │           ├── c1mim_0
    │   │           │   └── c1mim_0.xyz
    │   │           ├── nh3_1
    │   │           │   └── nh3_1.xyz
    │   │           ├── nh3_2
    │   │           │   └── nh3_2.xyz
    │   │           └── nh3_3
    │   │               └── nh3_3.xyz
    │   ├── ch
    │   │   └── ac
    │   │       ├── nowater
    │   │       │   ├── ch_ac_nowater.xyz
    │   │       │   └── frags
    │   │       │       ├── acetate_0
    │   │       │       │   └── acetate_0.xyz
    │   │       │       ├── acetate_1
    │   │       │       │   └── acetate_1.xyz
    │   │       │       ├── choline_2
    │   │       │       │   └── choline_2.xyz
    │   │       │       └── choline_3
    │   │       │           └── choline_3.xyz
    │   │       └── water
    │   │           ├── ch_ac_water.xyz
    │   │           └── frags
    │   │               ├── acetate_0
    │   │               │   └── acetate_0.xyz
    │   │               ├── acetate_1
    │   │               │   └── acetate_1.xyz
    │   │               ├── choline_2
    │   │               │   └── choline_2.xyz
    │   │               ├── choline_3
    │   │               │   └── choline_3.xyz
    │   │               └── water_4
    │   │                   └── water_4.xyz
    │   └── water
    │       ├── frags
    │       │   ├── water_0
    │       │   │   └── water_0.xyz
    │       │   └── water_1
    │       │       └── water_1.xyz
    │       └── water.xyz
    └── files
        ├── c1mim_nh3.xyz
        ├── ch_ac_nowater.xyz
        ├── ch_ac_water.xyz
        └── water.xyz"""

    # calc_dir = make_parent_dir()
    cwd = os.getcwd()
    for file in files:
        orig_dir = os.path.dirname(file)
        if not logfile_in_dir(orig_dir):
            new_dirs = make_dir_list(file)
            for idx, d in enumerate(new_dirs):
                # if dir is rerun, no need to move to lower dir
                change_to_subdir(d)
                if idx + 1 == len(new_dirs):  # when at maximum depth
                    copy_xyz(xyz_dir, file)
            os.chdir(cwd)


def xyz_is_rerun(file):
    return file == "rerun.xyz"


def logfile_in_dir(path="."):
    if path == "":
        path = "."
    return any("log" in f for f in os.listdir(path))


def make_job_files(base_dir, chem_package, settings):
    # find all xyz files in subdir to work on
    files = glob.glob("**/*xyz", recursive=True)
    for file in files:
        path, f = os.path.split(file)
        if path != "":  # or xyz_is_rerun(f):
            subdir = base_dir + "/" + path
            os.chdir(subdir)
            try:
                # print(f"Creating inputs for {file}...")
                # job = job_type(chem_package, f, settings)
                # job.create()
                job_type(chem_package, f, settings)
            except AttributeError as e:
                print(f">>> Error <<<")
                print(e)
            os.chdir(base_dir)


def make_job_subdirs(base_dir):
    """
    Look for input files in any subdirectory of ``calcs``, then creates a directory of that type
    (i.e. opt, spec, freq), then moves the inp and job into that folder.
    """
    parent = os.getcwd()
    for path, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".inp") or file.endswith(".job"):
                os.chdir(path)
                file_type = file[:-4]  # opt, spec, freq...
                os.mkdir(file_type)
                os.system(f"mv {file_type}.inp {file_type}.job {file_type}/")
                os.chdir(parent)


def xyz_to_tree(settings):
    """
    Takes a directory containing xyz files and creates a directory tree based on the filenames of
    the xyz files present. Uses underscores as delimiters for new subdirectories i.e. every time an
    underscore is seen, a new subdirectory is created.

    This function looks for a directory called ``files``, containing the xyz files, and outputs into
    ``calcs``. Ideally, call the function from the parent directory of ``files``. Desired settings
    should be created and then passed into the function.

    Another desirable feature is to separate optimisations from single point calculations. As a
    result, when the function runs, it creates jobs in a new subdirectory of the molecule
    directory.

    >>> s = Settings()
    >>> s.input.basis.gbasis = 'ccd' # gamess input (this is actually the default setting)
    >>> xyz_to_tree(s)

    Gives the following directory structure:
    .
    ├── calcs
    │   └── c1mim
    │       └── nh3
    │           ├── c1mim_nh3.xyz
    │           ├── frags
    │           │   ├── c1mim_0
    │           │   │   ├── c1mim_0.xyz
    │           │   │   └── opt
    │           │   │       ├── opt.inp
    │           │   │       └── opt.job
    │           │   ├── nh3_1
    │           │   │   ├── nh3_1.xyz
    │           │   │   └── opt
    │           │   │       ├── opt.inp
    │           │   │       └── opt.job
    │           │   ├── nh3_2
    │           │   │   ├── nh3_2.xyz
    │           │   │   └── opt
    │           │   │       ├── opt.inp
    │           │   │       └── opt.job
    │           │   └── nh3_3
    │           │       ├── nh3_3.xyz
    │           │       └── opt
    │           │           ├── opt.inp
    │           │           └── opt.job
    │           └── opt
    │               ├── opt.inp
    │               └── opt.job
    └── files
        └── c1mim_nh3.xyz

    Note: If a directory named ``calcs`` is already present, nonsensical results will be returned- any
    directory containing an xyz file will be acted upon. To run smoothly, remove or rename an existing
    ``calcs`` directory.
    """
    package = ask_package()
    xyz_directory = os.getcwd()
    files = get_xyz()
    files = [f for f in files if not logfile_in_dir(os.path.dirname(f))]
    make_tree_and_copy(xyz_directory, files)
    make_job_files(xyz_directory, package, settings)  # xyz directory is base dir
    os.chdir(xyz_directory)
