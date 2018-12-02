Python OOP project, implementing GAMESS, PSI4 and LAMMPS file production and
submission to SLURM and PBS scheduling systems. An important mention is that
this project focuses on implementing the Fragment Molecular Orbital approach to
quantum chemical calculations using GAMESS.

Objectives include:

- Automatic input and job file creation
- Scraping log files for relevant results; energies, coordinates, vibrations
- Logging of relevant results to a database for easy recall from any machine
  (credentials stored in ...)
- Reporting of errors resulting in failed calculations
- Automatic analysis of results (eventually, using dplython and plotnine to
  replicate the ease of  use found in dplyr and ggplot2 for R)

TODO: Add path of CLI to user's path when setup.py is run:
- export PATH=$PWD:$PATH  
