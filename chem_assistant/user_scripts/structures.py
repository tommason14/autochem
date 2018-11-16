def get_structures(basedir):
    """Returns a copy of the directory tree without the runtype folders for each chemical system- no need to know if the file was ran as an optimisation or hessian etc... the files with a chemical name are initial coordinates, and equil.xyz are equilibrium coordinates found after geometry optimisation."""

    # shutil.copytree
    parent = os.getcwd()
        for path, dirs, files in os.walk(base_dir):
            for file in files:
                pass
            
            


