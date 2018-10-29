class Input:
    """Base class for any input file for a computational chemistry calculation- ab initio or molecular dynamics.

    Instances of this class have the following attributes:
    * ``path`` -- file path of the desired input file
    * ``run`` -- type of calculation required
    * ``using`` -- coordinates of chemical system, in xyz format

    """
    def __init__(self, path, run, using = None):
        self.path = path
        self.run = run
        if using is not None:
            self.coords = read_xyz(using)
            #list of Atom objects, more useful than list of coordinates

    # def coords(self):
    #     """Reads in coordinates from an xyz file. Note this will be extended to work for mol and pdb file types in future"""
    #     atom =
    #     self.coords = []
    #     with open(self.coord_file, "r") as cf:
    #     for line in cf.readlines(): #memory efficiency
    #         if re.search(atom, line):
    #             coords.append(line)
    #     return self.coords


    def read_xyz(self, using):
        """Reads in coordinates from an xyz file"""
        coords = []
        with open(using, "r") as f:
            for coord in f.readlines()[2:]:
                line = coord.split()
                for val in PT.ptable.values():
                    if line[0] == val[0]:
                        coords.append(Atom(line[0], coords = tuple([float(i) for i in line[1:4]])))
        return coords

class Ab_initio(Input):
    """Base class for any input file for an ab initio calculation"""
    def __init__(self, path, run, using, basis, fmo_type = None, c_os = None):
        super.__init__(path, run, using)
        self.basis = basis
        if fmo_type is not None:
            self.fmo_type = fmo_type
        if c_os is not None:
            self.c_os = c_os


class Gamess_input(Ab_initio):
    """Class representing an input file for the GAMESS quantum chemistry package.

    Instances of this class will have the following attributes:

    * ``path`` -- file path to the generated input file
    * ``run`` -- type of calculation required. Can choose from the following options:
        * optimize
        * single point
        * hessian

        The default behaviour is to output files in accordance with Fragment Molecular Orbital theory. If no FMO data is needed, pass the option along with 'no fmo'
            >>> Gamess_input(..., run = "hessian no fmo", ...)
    * ``using`` -- file path pointing to the coordinates of the system. Only xyz files are supported currently, with plans to implement pdb and mol support in future.
    * ``basis`` -- basis set of the ab initio calculation.

    Usage:
        >>> Gamess_input(filename, run = "optimize", using = "file.xyz")
        >>> Gamess_input(filename, run = 'energy', using = "file.xyz", fmo_type=3)
    """

    # TODO: modify this class to work with GAMESS calcs for molecular dynamics. Make more options for run, set basis, fmo type, and c_os equal to none initially, and set as required for ab initio calcs

    def __init__(self, path, run, using, basis, fmo_type, c_os):
        super.__init__(path, run, using, basis, fmo_type, c_os)

        # avoid the need to call basis and c_os if just a normal run
        if c_os is None:
            if "opt" in self.run.lower():
                self.basis = 'CCD'
                self.c_os = 1.752
            elif "energy" in self.run_lower():
                self.basis = 'CCT'
                self.c_os = 1.64

    def header(self):
        return gamess_header[run]


    gamess_header = {}

    gamess_header['optimize'] = f"""\
 $SYSTEM MWORDS=500 MEMDDI=0 $END
 $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE MAXIT=200 ISPHER=1 $END
 $GDDI NGROUP=8 $END
 $STATPT NSTEP=500 $END
 $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
 $BASIS GBASIS={self.basis} $END
 $FMO
    NFRAG=8 NBODY=2
    MPLEVL(1)=2
    INDAT(1)=0,1,-25,
             0,26,-32,
             0,33,-39,
             0,40,-64,
             0,65,-67,
             0,68,-70,
             0,71,-73,
             0,74,-76,
             0
    ICHARG(1)=1,-1,-1,1,0,0,0,0
    RESPAP=0 RESPPC=-1 RESDIM=100 RCORSD=100
 $END
 $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
 $DATA
"""

    gamess_header['optimize no fmo'] = f"""\
 $SYSTEM MWORDS=600 $END
 $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE MAXIT=200 ISPHER=1 ICHARG=0 $END
 $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
 $BASIS GBASIS={self.basis} $END
 $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
 $DATA
"""


    gamess_header['single point'] = f"""\
 $SYSTEM MWORDS=500 MEMDDI=0 $END
 $CONTRL SCFTYP=RHF RUNTYP=ENERGY MAXIT=200 ISPHER=1 $END
 $GDDI NGROUP=4 $END
 $STATPT NSTEP=400 $END
 $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
 $BASIS GBASIS={self.basis} $END
 $FMO
    NFRAG=4 NBODY=1
    MPLEVL(1)=2
    INDAT(1)=0,1,-21,
             0,22,-28,
             0
    ICHARG(1)=1,-1
    RESPAP=0 RESPPC=-1 RESDIM=4 RCORSD=4
 $END
 $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
 $DATA
"""

    gamess_header['single point no fmo'] = f"""\
 $SYSTEM MWORDS=600 $END
 $CONTRL SCFTYP=RHF RUNTYP=ENERGY MAXIT=200 ISPHER=1 ICHARG=0 $END
 $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
 $BASIS GBASIS={self.basis} $END
 $MP2 CODE=IMS SCSPT=SCS SCSOPO={self.c_os} SCSPAR=0.0 $END
 $DATA
"""

    gamess_header['hessian'] = ""
    gamess_header['hessian no fmo'] = f"""\
 $SYSTEM MWORDS=1000 MEMDDI=0 $END
 $CONTRL SCFTYP=RHF MPLEVL=2 RUNTYP=HESSIAN MAXIT=200 ISPHER=1 $END
 $STATPT NSTEP=500 $END
 $BASIS GBASIS=CCD $END
 $SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. $END
 $MP2 CODE=IMS SCSPT=SCS SCSOPO=1.752 SCSPAR=0.0 $END
 $DATA
"""
