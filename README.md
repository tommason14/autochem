Python OOP project, implementing GAMESS, GAUSSIAN, PSI4 and ORCA file production and
submission to SLURM and PBS scheduling systems. An important point is that
this project focuses on implementing the Fragment Molecular Orbital approach to
quantum chemical calculations using GAMESS, and input files are generated for
SRS-MP2 jobs by default in GAMESS.

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
