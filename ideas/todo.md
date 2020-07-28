# Results
 
## Properties

- basis
- scf_data
- mp2_data
- homo
- lumo
- homo_lumo_gap

## Functions

- get_data()


# core.thermo file is not needed. Refactor so that the Results classes extract
# that data and the scripts.thermochemistry() function can write the freq.out
# file needed by the fortran code.
