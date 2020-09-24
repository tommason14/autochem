# autochem

# Overview

Python OOP project, implementing GAMESS, GAUSSIAN, PSI4 and ORCA file production and
submission to SLURM and PBS scheduling systems. GAMESS implementations focus
heavily on the Fragment Molecular Orbital approach to
quantum chemical calculations, and input files are generated for
[SRS-MP2](https://aip.scitation.org/doi/10.1063/1.4975326) jobs by default in GAMESS.

Use this code for:

- Automatic input and job file creation
- Scraping log files for relevant results: 
  - energies
  - geometries 
    - intermolecular hydrogen bond lengths
  - vibrations
  - fluorescence data
  - geodesic charges
  - homo-lumo gaps
- Automatic analysis of results
  - interaction energies
    - purely ionic systems
    - mixed ionic/neutral species
  - calculates free energies

Also exportable as a python package, and can be easily extended. See the
[examples](https://github.com/tommason14/monash_automation/tree/master/examples) for more details.

Currently set up for the following supercomputers:
- Gadi (Canberra)
- Raijin (Canberra, now decommissioned)
- Magnus (Perth)
- Monarch (Melbourne)
- Massive (Melbourne)
- Stampede2 (Texas)

# Installation

Pip install available: `pip install autochem`
Beware, this version may have bugs and is only updated semi-regularly.

Recommended installation (settings files are easier to access): 
```sh
git clone https://github.com/tommason14/autochem
cd autochem
sh install.sh
```
  
# Example Usage

# Creating input files for computational calculations

Important: When creating job files, all data for a job can be modified through the use of a
`Settings` object.

Users have various options when choosing how to create input and job files:

1. Create a python script including the desired input file and job parameters,
   and run the job using a python class. For example:

```python
from autochem import GaussJob, Settings
from glob import glob

xyz = glob('*xyz')[0] # find the xyz file to run

sett=Settings()
# define the input parameters
sett.input.method='wb97xd'
sett.input.basis='aug-cc-pVDZ'
sett.input.opt=True
sett.input.td='nstates=10,root=9'
sett.input.scrf='smd,solvent=n,n-DiMethylFormamide'
# define job scheduler parameters
sett.meta.mem='160gb'
sett.meta.ncpus=48
sett.meta.nodemem='192GB' # gaussian only
sett.meta.time='24:00:00'
sett.meta.jobfs='400gb' # gadi only

# create the input and job files
GaussJob(xyz, settings=sett)
```
Using this method, options such as `frags_in_subdir=True` can be used, or
`filename` can be used to give the name of input/job files. Otherwise, jobs are
named according to their run types (optimisations are named `opt`, single points
are named `spec` and frequency calculations are named `hess`).

2. Store the parameters of the `Settings` object in a python file, then run from
   the command line:

```python
# settings.py
sett.input.method='wb97xd'
sett.input.basis='aug-cc-pVDZ'
sett.input.opt=True
sett.input.td='nstates=10,root=9'
sett.input.scrf='smd,solvent=n,n-DiMethylFormamide'

sett.meta.mem='160gb'
sett.meta.ncpus=48
sett.meta.nodemem='192GB'
sett.meta.time='24:00:00'
sett.meta.jobfs='400gb'
```
then run from the command line using `autochem -d -s settings.py`, which
takes in every xyz file in the current directory and creates jobs using
parameters from the `Settings` object in `settings.py`. 
Note: using this method, the `Settings` object must be called `sett`.

Example settings can be found in the [settings_files](https://github.com/tommason14/monash_automation/tree/master/settings_files)
directory.

Methods of defining input and job parameters for each software are different, and
outlined below.

## GAMESS-US

Using a dummy system, `h.xyz`, we can see the default options like so:

```python
from autochem import GamessJob

gamess = GamessJob('h.xyz')

print(gamess.input)
```

which produces:

```
basis:    
      gbasis:    ccd
contrl:    
       icharg:    0
       ispher:    1
       maxit:    200
       mplevl:    2
       runtyp:    optimize
       scftyp:    rhf
mp2:    
    code:    ims
    scsopo:    1.752
    scspar:    0.0
    scspt:    scs
scf:    
    diis:    .true.
    dirscf:    .true.
    fdiff:    .false.
statpt:    
       nstep:    500
system:    
       memddi:    0
       mwords:    500
```

Note that the default settings produce an SRS-MP2 optimisation using a cc-pVDZ
basis set. Here is the resulting input file:

```
 $SYSTEM MEMDDI=0 MWORDS=500 $END
 $CONTRL ICHARG=0 ISPHER=1 MAXIT=200 MPLEVL=2 
  RUNTYP=OPTIMIZE SCFTYP=RHF $END
 $STATPT NSTEP=500 $END
 $SCF DIIS=.TRUE. DIRSCF=.TRUE. FDIFF=.FALSE. $END
 $BASIS GBASIS=CCD $END
 $MP2 CODE=IMS SCSOPO=1.752 SCSPAR=0.0 SCSPT=SCS $END
 $DATA
h
C1
 H       1.0   1.00000    2.00000    3.00000
 $END
```

To change the parameters of the file, we use a `Settings` object.

To run an open shell single point calculation:

```python
from autochem import Settings, GamessJob

sett = Settings()
sett.input.contrl.runtyp='energy'
sett.input.contrl.scftyp='ROHF'

gamess = GamessJob('h.xyz', settings=sett)
```

Giving us this file:
```
 $SYSTEM MEMDDI=0 MWORDS=500 $END
 $CONTRL ICHARG=0 ISPHER=1 MAXIT=200 MPLEVL=2 
  RUNTYP=ENERGY SCFTYP=ROHF $END
 $SCF DIIS=.TRUE. DIRSCF=.TRUE. FDIFF=.FALSE. $END
 $BASIS GBASIS=CCD $END
 $MP2 CODE=IMS SCSOPO=1.752 SCSPAR=0.0 SCSPT=SCS $END
 $DATA
h
C1
 H       1.0   1.00000    2.00000    3.00000
 $END
```

### Writing your own `Settings` objects

We are modifying the input file here, so use
`sett.input`. Then as GAMESS uses a `$FLAG KEYWORD=VALUE` syntax, that is
implemented here as `sett.input.flag.keyword = value`.

To generate a GEODESIC partial charge calculation at the Hartree-Fock level,
we have to include some parameters to the `$ELPOT` (electric potential)
section, and remove the MP2 default settings:

```python
from autochem import Settings

sett=Settings()
sett.input.basis.gbasis='cct'
sett.input.contrl.runtyp='energy'
sett.input.elpot.iepot=1
sett.input.elpot.where='pdc'
sett.input.pdc.ptsel='geodesic'

# rm defaults
sett.input.contrl.mplevl=None
sett.input.mp2=None
```

This will then generate the following file:

```
 $SYSTEM MEMDDI=0 MWORDS=500 $END
 $CONTRL ICHARG=0 ISPHER=1 MAXIT=200 RUNTYP=ENERGY SCFTYP=RHF $END
 $ELPOT IEPOT=1 WHERE=PDC $END
 $PDC PTSEL=GEODESIC $END
 $SCF DIIS=.TRUE. DIRSCF=.TRUE. FDIFF=.FALSE. $END
 $BASIS GBASIS=CCT $END
 $DATA
title
C1
 coords...
 $END
```

To run DFT calculations, remove the MP2 defaults and then add the appropriate
settings:

```python
from autochem import Settings

sett=Settings()
sett.input.mp2=None
sett.input.contrl.mplevl=None
sett.input.contrl.dfttyp='m06-2x'
sett.input.dft.method='grid'
```

### Additional options

When running FMO calculations, you may wish to group certain molecules together.
This can be acheived using `sett.grouped='sodium-bf4`, to group together sodium
and tetrafluoroborate ions into one fragment.

Alternatively, large molecules could be too large to run as one fragment. In
that case, use `sett.bonds_to_split=[(28,29), (40,41)]`, and pass in a nested
list of atoms that form the bonds that should be broken. (Currently
experimental and may not work as desired.)

FMO jobs are run by using the `GamessJob(..., fmo=True)` option. If running
using the command line (`autochem -d`), FMO jobs can also be chosen.

### Job information

Information for the SLURM/PBS schedulers are given as `sett.meta.option=choice`.
For GAMESS, options include `ncpus`, `mem`, `partition`, `time` and `jobfs` (PBS
only).
If you are submitting a job on monarch or massive, a `sett.meta.nodes`
parameter is also used. For example:

```
sett.supercomp='mas'
sett.meta.ncpus=60
sett.meta.mem=120
sett.meta.nodes=2
```

produces:

```
#!/bin/env bash
#SBATCH --job-name=opt
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --ntasks=60
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=30
...
```

To give information on job information for each fragment when choosing
`GamessJob(..., frags_in_subdir=True)`, use `sett.frag.meta.option=choice`.

## GAUSSIAN

GAUSSIAN commands are defined in groups with a `keyword` or `keyword=value`
syntax, such as `opt`, `opt=loose`, or `opt=(calcfc,noeigentest,ts)`.

This programs achieves the desired output by looking at the parameters set
inside the `Settings` object and deciding which output to choose accordingly.

The following `Settings` file:

```python
from autochem import Settings

sett=Settings()
sett.input.opt='ts,noeigentest,calcfc'
sett.input.freq=True
sett.input.scrf='smd,solvent=water'
sett.input.charge=0
sett.input.mult=2
```

Produces the following input file:

```
%chk=opt-freq.chk
%mem=160gb
%nproc=48

#P M062X/cc-pVDZ opt=(ts,noeigentest,calcfc) freq int=(grid=ultrafine) scrf=(smd,solvent=water)

rerun

0 2
coords...
```

### Job information

Information for the SLURM/PBS schedulers are given as `sett.meta.option=choice`.
For GAUSSIAN, options include `ncpus`, `mem` (%mem=... in the input file),
`nodemem` (defining memory for the scheduler), `partition`, `time` and `jobfs`
(PBS only).

To give information on job information for each fragment when choosing
`GaussJob(..., frags_in_subdir=True)`, use `sett.frag.meta.option=choice`.

## PSI4

PSI4 options are given as:
- `self.input.molecule` for parameters inside the molecule section:
```
molecule {
charge mult
...
units ...
symmetry ...
}
```
- `sett.input.globals` for parameters inside the globals section:
```
set globals {
    basis ...
    scf_type ...
}
```
The default settings are as follows:
```python
from autochem import PsiJob
psi = PsiJob('file.xyz')
print(psi.input)
```

```
charge:    0
globals:
        S_ORTHOGONALIZATION:    canonical
        basis:    cc-pVTZ
        freeze_core:    True
        guess:    sad
        scf_type:    DF
memory:    30 gb
molecule:
         charge:    0
         multiplicity:    1
         symmetry:    c1
         units:    angstrom
mult:    1
run:
    energy:    mp2
unbound:
```

To change the run type, use `sett.input.run = {'optimize': 'scf'}` to produce
`optimize('scf')`.
For additional run types, use `sett.input.run.additional = {'dertype': 'energy'}` 
to produce `optimize('scf', dertype='energy')`. 

A possible `Settings` object may look like:
```python
from autochem import Settings

sett = Settings()
sett.input.memory='60gb'
sett.input.globals.basis='aug-cc-pVTZ'
sett.input.charge=2
sett.input.mult=3
sett.input.run = {'optimize': 'scf'}
sett.input.run.additional = {'dertype': 'energy'}
```
which produces:
```
# PSI4 Calc: xyz

memory 60gb

molecule complex {
2 3
 C       -3.94373    0.13147    0.12554
 C       -2.80294   -0.70993    0.69713
 ...
units angstrom
no_reorient
symmetry c1
}

set globals {
    S_ORTHOGONALIZATION canonical
    basis aug-cc-pVTZ
    freeze_core True
    guess sad
    scf_type DF
}
optimize('scf', dertype='energy')
```

### Counterpoise correction

To produce counterpoise corrected jobs, use `PsiJob(..., cp=True)`.

### Job information

Information for the SLURM/PBS schedulers are given as `sett.meta.option=choice`.
For PSI4, options include `ncpus`, `mem`, `time`, `partition`, and `jobfs` (PBS
only).

To give information on job information for each fragment when choosing
`PsiJob(..., frags_in_subdir=True)`, use `sett.frag.meta.option=choice`.

## ORCA

ORCA jobs rely on settings given as `!wB97X-D3 aug-cc-pVDZ RIJCOSX`, 
along with more detailed commands given as
```
%pal
  nprocs 48
end
```
to produce a job that will run on 48 cpus.

This is acheived by the following `Settings` file:

```python
from autochem import Settings

sett=Settings()
sett.input.run='wB97X-D3 aug-cc-pVTZ RIJCOSX'
sett.input.meta.tddft="""\
  nroots 10"""
sett.input.meta.cpcm="""\
  SMD true
  SMDSolvent "DMSO" """
```

Which produces the following input file:

```
!wB97X-D3 aug-cc-pVDZ RIJCOSX CPCM

%pal
 nprocs 48
end

%cpcm
  SMD true
  SMDSolvent "DMSO"
end

%tddft
  nroots 10
  tda false
end

*xyzfile ...
```

### Job information

Information for the SLURM/PBS schedulers are given as `sett.meta.option=choice`.
For ORCA, options include `ncpus`, `mem`, `time`, `partition`, and `jobfs` (PBS
only).

To give information on job information for each fragment when choosing
`OrcaJob(..., frags_in_subdir=True)`, use `sett.frag.meta.option=choice`.

## Inputs for all molecules in an xyz file

All calculations have the option of producing files for each molecule in the
system. To do this, pass in a `frags_in_subdir` option. For example,
`GamessJob('file.xyz', frags_in_subdir=True, settings=sett)`.

When running from the command line with `autochem -d`, input files can 
also be created for all molecules in the system.

In `Settings` files, to specify job parameters type `sett.frag.meta...`, whereas
job parameters for the entire complex are given as `self.meta...`

Note that an ionic cluster will also be produced i.e. the original xyz file 
with all neutral molecules removed.

## Information for job schedulers

The program is designed to run on remote supercomputers, and uses the `hostname`
of the user's account to decide which job file should be created. Alternatively,
users can pass in a `sett.supercomp='cluster'` option to force the program to create job
scripts for a desired cluster. Options include 'gadi', 'rjn' or 'raijin', 'mgs'
or 'magnus', 'mon' or 'monarch', 'm3', 'mas' or 'massive', and 'stm' or
'stampede'.
Any of those options will work.
For GAMESS jobs, separate job scripts are provided to use either the standard
versions of GAMESS, or the modified version that implements the option to choose
MP2 spin parameters.

Submission of jobs can be handled by a shell script such as:

```bash
#!/bin/bash

[[ $HOSTNAME =~ gadi ]] && submit="qsub" || submit="sbatch"

cwd=$(pwd)
for f in $(find . -name "*job")
do
    cd "$(dirname $f)"
    logs_in_dir=$(ls *log 2>/dev/null | wc -l)
    [[ $logs_in_dir -eq 0 ]] && $submit "$(basename $f)"
    cd $cwd
done
```

# Obtaining results

This program allows users to extract data using a command line interface.
The table below highlights the desired data and the command required:

Data  | Command
:---: | :---:
Energies | `autochem -r`
Partial charges | `autochem --charges`
Fluorescence data | `autochem --fluorescence`
HOMO-LUMO data | `autochem --homo-lumo`
Hydrogen bond data | `autochem -b`
Thermochemistry data | `autochem -t [temp in K] -m [multiplicity]`
Frequencies | `autochem --freqs-to-csv`

Every command shown above produces a csv file. These commands also 
allow you to give a filename with the `-o` parameter. For example, 
`autochem --homo-lumo -o homo_lumo.csv` to save the data into
`homo_lumo.csv`.

Outputs can also be limited to files that contain
a certain string in their path. To do this, use the `-l` flag. For example,
`autochem -rl 'spec'` to return only single point energies.

In addition, other information can be found:
- to look for equilibrated coordinates, use `autochem -e`
  - this will create either `spec` or `rerun` directories as subdirectories of
    each log file, depending on whether the optimisations have completed.
- to print out fragments of each xyz file in a directory, use `autochem -p`.
  If you wish for fragments to be grouped, use 
  `autochem -p -g 'lithium-sacchrinate'`. To give a more verbose output,
  showing the atom numbers of each fragment, use the `-v` flag.
- Interaction energies can be calculated using the output of `autochem -r`,
  by using the `-c` flag.
  By default, this assumes that you want to calculate the interaction of each
  species in the cluster. i.e. cluster - sum(all_molecules).
  If you wish to calculate the interaction of neutral species i.e. cluster -
  ionic - sum(neutral), pass in a `-with-ionic` flag. For example, 
  `autochem -c [results.csv] --with-ionic`. (Experimental, use with caution.)
- Gibbs free energies can also be calculated, by taking in the thermochemical
  data of `autochem -t [temp in K] -m [multiplicity]`, the interaction
  energies of `autochem -c [results.csv] [--with-ionic]`. (Experimental, use
  with caution.)
- Boltzmann-weighted interaction energies can be calculated using the output of
  `autochem -c [results.csv] [--with-ionic]`. A grouping parameter is
  required here, and is written as though the csv file is described as `df`, a
  `pandas.DataFrame` object. For example, 
  `autochem -w data.csv --group df['Config'].str.split('-').str[:-1].str.join('-')`. 
  (Experimental, use with caution)
  
# Adding additional molecules to the database

Add molecules to the `~/.config/autochem/molecules.txt` file, 
using the format of:
```
name=dihydrogen_citrate
charge=-1
multiplicity=1
atoms=C,C,C,H,H,O,O,C,O,C,C,H,H,O,H,O,H,H,O,O
```

Make sure that names do not contain spaces.
