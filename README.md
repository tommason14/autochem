Python OOP project, implementing GAMESS, GAUSSIAN, PSI4 and ORCA file production and
submission to SLURM and PBS scheduling systems. An important point is that
this project focuses on implementing the Fragment Molecular Orbital approach to
quantum chemical calculations using GAMESS, and input files are generated for
SRS-MP2 jobs by default in GAMESS.

Use this code for:

- Automatic input and job file creation
- Scraping log files for relevant results: energies, coordinates, vibrations,
  fluorescence data, geodesic charges
- Automatic analysis of results
  - Interaction and free energies
