import os

__all__ = ['make_files_from_meta']

def make_files_from_meta(base_dir, filename = None):
    """This function looks for ``meta.py`` in any subdirectory below the directory passed as an argument. ``meta.py`` is then run as a python file, with the idea of including data for a job. As long as there is an xyz file in the same directory as ``meta.py``,  the desired result is passed. 
    
    Usage:
        Files in directory:
            - ch_ac.xyz
            - meta.py
        
        meta.py:

            from autochem import GamessJob
            GamessJob(using='ch_ac.xyz', frags_in_subdir = True)

    Resulting directory structure:
    
        dir
        ├── ch_ac.xyz
        ├── frags
        │   ├── acetate_0
        │   │   ├── acetate_0.xyz
        │   │   └── opt
        │   │       ├── opt.inp
        │   │       └── opt.job
        │   ├── acetate_1
        │   │   ├── acetate_1.xyz
        │   │   └── opt
        │   │       ├── opt.inp
        │   │       └── opt.job
        │   ├── choline_2
        │   │   ├── choline_2.xyz
        │   │   └── opt
        │   │       ├── opt.inp
        │   │       └── opt.job
        │   ├── choline_3
        │   │   ├── choline_3.xyz
        │   │   └── opt
        │   │       ├── opt.inp
        │   │       └── opt.job
        │   └── water_4
        │       ├── opt
        │       │   ├── opt.inp
        │       │   └── opt.job
        │       └── water_4.xyz
        ├── meta.py
        └── opt
            ├── opt.inp
            └── opt.job
    """
    if filename is not None:
        parent = os.getcwd()
        for path, dirs, files in os.walk(base_dir):
            for file in files:
                if file == 'meta.py' and os.path.exists(os.path.join(path, 'equil.xyz')):
                    print(path)
                    os.chdir(path)
                    if not any(file.endswith('.xyz') for file in os.listdir('.')):
                        raise TypeError(f'Meta file requires an xyz file in the same directory. Check {os.getcwd()}')
                    else:
                        print(os.getcwd())
                    os.system('python3 meta.py')
                    os.chdir(parent)

    else:
        parent = os.getcwd()
        for path, dirs, files in os.walk(base_dir):
            for file in files:
                if file == 'meta.py':
                    os.chdir(path)
                    if not any(file.endswith('.xyz') for file in os.listdir('.')):
                        raise TypeError(f'Meta file requires an xyz file in the same directory. Check {os.getcwd()}')
                    else:
                        print(os.getcwd())
                    os.system('python3 meta.py')
                    os.chdir(parent)

    
        
            



