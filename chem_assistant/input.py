from atom import Atom
from molecule import Molecule
from settings import *

class Input:
    """Base class for any input file for a computational chemistry calculation- ab initio or molecular dynamics.

    Instances of this class have the following attributes:
    * ``path`` -- file path of the desired input file
    * ``run`` -- type of calculation required
    * ``using`` -- coordinates of chemical system, in xyz format

    """
    def __init__(self, using = None):
        if using is not None:
            self.mol = Molecule(using)
            #self.mol.separate() # in the "user side" file
            #self.mol.gamess_format() # in the "user side" file
            # for frag in self.mol.fragments: Input_type(using=f'fragments/{frag}') -> no fmo though

class Gamess_input(Input):
    
    def __init__(self, using = None, fmo = False):
        super().__init__(using)
        if fmo:
            self.input = read_template('gamess_fmo.json')
            self.mol.separate()
            self.ind, self.ch = self.mol.gamess_format
        else:
            self.input = read_template('gamess.json')


# class Ab_initio(Input):
#     """Base class for any input file for an ab initio calculation"""
#     def __init__(self, path, run, using, basis, fmo_type = None, c_os = None):
#         super.__init__(path, run, using)
#         self.basis = basis
#         if fmo_type is not None:
#             self.fmo_type = fmo_type
#         if c_os is not None:
#             self.c_os = c_os


# class Gamess_input(Ab_initio):
#     """Class representing an input file for the GAMESS quantum chemistry package.

#     Instances of this class will have the following attributes:

#     * ``path`` -- file path to the generated input file
#     * ``run`` -- type of calculation required. Can choose from the following options:
#         * optimize
#         * single point
#         * hessian

#         The default behaviour is to output files in accordance with Fragment Molecular Orbital theory. If no FMO data is needed, pass the option along with 'no fmo'
#             >>> Gamess_input(..., run = "hessian no fmo", ...)
#     * ``using`` -- file path pointing to the coordinates of the system. Only xyz files are supported currently, with plans to implement pdb and mol support in future.
#     * ``basis`` -- basis set of the ab initio calculation.

#     Usage:
#         >>> Gamess_input(filename, run = "optimize", using = "file.xyz")
#         >>> Gamess_input(filename, run = 'energy', using = "file.xyz", fmo_type=3)
#     """

#     # TODO: modify this class to work with GAMESS calcs for molecular dynamics. Make more options for run, set basis, fmo type, and c_os equal to none initially, and set as required for ab initio calcs

#     def __init__(self, path, run, using, basis, fmo_type, c_os):
#         super.__init__(path, run, using, basis, fmo_type, c_os)

#         # avoid the need to call basis and c_os if just a normal run
#         if c_os is None:
#             if "opt" in self.run.lower():
#                 self.basis = 'CCD'
#                 self.c_os = 1.752
#             elif "energy" in self.run_lower():
#                 self.basis = 'CCT'
#                 self.c_os = 1.64

#     def header(self):

#         run_type = ''
#         if 'opt' in self.run:
#             run_type = 'optimize'
#         elif 'single' in self.run or 'sp' in self.run:
#             run_type = 'single point'
#         elif 'hess' in self.run or 'freq' in self.run:
#             run_type = 'hessian'

#         # if no self.fmo_type defined, use no fmo versions        
#         if not hasattr(self, fmo_type):
#             run_type += ' no fmo'

#         return gamess_header[run_type]



#     gamess_header = {}

#     gamess_header['optimize'] = f"""\
#  $SYSTEM MWORDS=500 MEMDDI=0 $END
#  $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE MAXIT=200 ISPHER=1 ICHARG={self.complex['charge']} $END
#  $GDDI NGROUP={len(self.mol.fragments)} $END
#  $STATPT NSTEP=500 $END
#  $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
#  $BASIS GBASIS={self.basis} $END
#  $FMO
#     NFRAG={len(self.mol.fragments)} NBODY=2
#     MPLEVL(1)=2
#     INDAT(1)={self.mol.indat}
#     ICHARG(1)={self.mol.charg}
#     RESPAP=0 RESPPC=-1 RESDIM=100 RCORSD=100
#  $END
#  $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
#  $DATA
# """

#     gamess_header['optimize no fmo'] = f"""\
#  $SYSTEM MWORDS=500 $END
#  $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE MAXIT=200 ISPHER=1 ICHARG=charge $END
#  $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
#  $BASIS GBASIS={self.basis} $END
#  $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
#  $DATA
# """


#     gamess_header['single point'] = f"""\
#  $SYSTEM MWORDS=500 MEMDDI=0 $END
#  $CONTRL SCFTYP=RHF RUNTYP=ENERGY MAXIT=200 ISPHER=1 ICHARG={self.complex['charge']}$END
#  $GDDI NGROUP={len(self.mol.fragments)} $END
#  $STATPT NSTEP=400 $END
#  $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
#  $BASIS GBASIS={self.basis} $END
#  $FMO
#     NFRAG={len(self.mol.fragments)} NBODY=1
#     MPLEVL(1)=2
#     INDAT(1)={self.mol.indat}
#     ICHARG(1)={self.mol.charg}
#     RESPAP=0 RESPPC=-1 RESDIM=100 RCORSD=100
#  $END
#  $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
#  $DATA
# """

#     gamess_header['single point no fmo'] = f"""\
#  $SYSTEM MWORDS=500 $END
#  $CONTRL SCFTYP=RHF RUNTYP=ENERGY MAXIT=200 ISPHER=1 ICHARG=charge $END
#  $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
#  $BASIS GBASIS={self.basis} $END
#  $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
#  $DATA
# """

#     gamess_header['hessian'] = ""
#     gamess_header['hessian no fmo'] = f"""\
#  $SYSTEM MWORDS=500 MEMDDI=0 $END
#  $CONTRL SCFTYP=RHF MPLEVL=2 RUNTYP=HESSIAN MAXIT=200 ISPHER=1 ICHARG=charge $END
#  $STATPT NSTEP=500 $END
#  $BASIS GBASIS={self.basis} $END
#  $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
#  $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
#  $DATA
# """
#     # build gamess_input (gi)
#     gi = gamess_header[run]
#     gi += f"{self.run}--{using[:-4]}" #strip.xyz from coordinate file name #title
#     gi += 'C1\n'
#     for element in self.mol.formula(as_dict = True):    # N 7.0
#         gi += f" {element} {PT.get_atnum(element)}.0\n"  # S 16.0
#     gi += " $END\n"
#     for atom in self.coords:
#     gi += f" {atom.symbol:5s} {PT.get_atnum(atom.symbol)}.0 {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}\n"
#     gi += ' $END'



