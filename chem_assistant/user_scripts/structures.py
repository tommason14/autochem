import os

__all__ = ['copy_xyz_tree']

def get_structures(base_dir):
    """Returns a copy of the directory tree without the runtype folders for each chemical system- no need to know if the file was ran as an optimisation or hessian etc... the files with a chemical name are initial coordinates, and equil.xyz are equilibrium coordinates found after geometry optimisation."""

    file_lst = []
    for path, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.xyz'):
                file_lst.append(os.path.join(path, file))
    return file_lst

def change_to_subdir(subdirectory):
    new_dir = os.path.join(os.getcwd(), subdirectory)
    if os.path.isdir(new_dir):
        os.chdir(new_dir)
    else: 
        os.mkdir(new_dir)
        os.chdir(new_dir)

def make_dir_list(file): 
    new_dirs = []
    for part in file.split('/'):
        new_dirs.append(part)
    return new_dirs

def copy_tree(xyzlist, base_dir, new_dir):
    """
    Copies the directory structure from the base directory downwards, and then when the final
    subdirectory has been created, files are copied over
    """
    base = os.path.abspath(base_dir)
    new = os.path.join(base, new_dir)
    if not os.path.exists(new):
        os.mkdir(new)
    os.chdir(new)
    parent = os.getcwd()
    for file in xyzlist:
        filepath, filename = os.path.split(file)
        new_dirs = make_dir_list(filepath)
        for idx, d in enumerate(new_dirs):
            change_to_subdir(d)
            if idx + 1 == len(new_dirs):
                # at bottom of dir tree, copy xyz over
                orig = os.path.join(base, filepath, filename)
                new = os.path.join(os.getcwd(), filename)
                os.system(f'cp {orig} {new}')
        os.chdir(parent)
    

def copy_xyz_tree(base_dir, new_dir):
    """
    Copy all xyz files recursively to a new directory, maintaining the original 
    directory structure. Give the new directory as a relative path to the 
    base directory that you are copying from.
    """
    structures = get_structures(base_dir)
    copy_tree(structures, base_dir, new_dir)
