# Monash Automation

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
[examples](examples/) section for more details.

# Example Usage

When creating job files, all data for a job can be modified through the use of a
`Settings` object.

Using a dummy system, `h.xyz`, we can see the default options like so:

```{python}
from chem_assistant import GamessJob

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

```{python}
from chem_assistant import Settings, GamessJob

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

# Writing your own `Settings` objects

We are modifying the input file here, so use
`sett.input`. Then as GAMESS uses a `$FLAG KEYWORD=VALUE` syntax, that is
implemented here as `sett.input.flag.keyword = value`.