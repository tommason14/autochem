from errors import PTError

__all__ = ['PeriodicTable']

class PeriodicTable:
    """Helper class to allow for lookup of atomic properties. Can convert between symbol and atomic number"""
    ptable = {}
    #[symbol, mass, radius, connectors]
        #atomic weights from: http://www.ciaaw.org/atomic-weights.htm
    ptable[  0] = ['Xx',   0.00000, 0.00 ,  0]
    ptable[  1] = [ 'H',   1.00798, 0.30 ,  1]
    ptable[  2] = ['He',   4.00260, 0.99 ,  0]
    ptable[  3] = ['Li',   6.96750, 1.52 ,  8]
    ptable[  4] = ['Be',   9.01218, 1.12 ,  8]
    ptable[  5] = [ 'B',  10.81350, 0.88 ,  6]
    ptable[  6] = [ 'C',  12.01060, 0.77 ,  4]
    ptable[  7] = [ 'N',  14.00685, 0.70 ,  3]
    ptable[  8] = [ 'O',  15.99940, 0.66 ,  2]
    ptable[  9] = [ 'F',  18.99840, 0.64 ,  1]
    ptable[ 10] = ['Ne',  20.17970, 1.60 ,  0]
    ptable[ 11] = ['Na',  22.98977, 1.86 ,  8]
    ptable[ 12] = ['Mg',  24.30550, 1.60 ,  8]
    ptable[ 13] = ['Al',  26.98154, 1.43 ,  8]
    ptable[ 14] = ['Si',  28.08500, 1.17 ,  8]
    ptable[ 15] = [ 'P',  30.97376, 1.10 ,  8]
    ptable[ 16] = [ 'S',  32.06750, 1.04 ,  2]
    ptable[ 17] = ['Cl',  35.45150, 0.99 ,  1]
    ptable[ 18] = ['Ar',  39.94800, 1.92 ,  0]
    ptable[ 19] = [ 'K',  39.09830, 2.31 ,  8]
    ptable[ 20] = ['Ca',  40.07800, 1.97 ,  8]
    ptable[ 21] = ['Sc',  44.95591, 1.60 ,  8]
    ptable[ 22] = ['Ti',  47.86700, 1.46 ,  8]
    ptable[ 23] = [ 'V',  50.94150, 1.31 ,  8]
    ptable[ 24] = ['Cr',  51.99610, 1.25 ,  8]
    ptable[ 25] = ['Mn',  54.93804, 1.29 ,  8]
    ptable[ 26] = ['Fe',  55.84500, 1.26 ,  8]
    ptable[ 27] = ['Co',  58.93319, 1.25 ,  8]
    ptable[ 28] = ['Ni',  58.69340, 1.24 ,  8]
    ptable[ 29] = ['Cu',  63.54600, 1.28 ,  8]
    ptable[ 30] = ['Zn',  65.38000, 1.33 ,  8]
    ptable[ 31] = ['Ga',  69.72300, 1.41 ,  8]
    ptable[ 32] = ['Ge',  72.63000, 1.22 ,  8]
    ptable[ 33] = ['As',  74.92159, 1.21 ,  8]
    ptable[ 34] = ['Se',  78.97100, 1.17 ,  8]
    ptable[ 35] = ['Br',  79.90400, 1.14 ,  1]
    ptable[ 36] = ['Kr',  83.79800, 1.97 ,  0]
    ptable[ 37] = ['Rb',  85.46780, 2.44 ,  8]
    ptable[ 38] = ['Sr',  87.62000, 2.15 ,  8]
    ptable[ 39] = [ 'Y',  88.90584, 1.80 ,  8]
    ptable[ 40] = ['Zr',  91.22400, 1.57 ,  8]
    ptable[ 41] = ['Nb',  92.90637, 1.41 ,  8]
    ptable[ 42] = ['Mo',  95.95000, 1.36 ,  8]
    ptable[ 43] = ['Tc',  98.00000, 1.35 ,  8]
    ptable[ 44] = ['Ru', 101.07000, 1.33 ,  8]
    ptable[ 45] = ['Rh', 102.90550, 1.34 ,  8]
    ptable[ 46] = ['Pd', 106.42000, 1.38 ,  8]
    ptable[ 47] = ['Ag', 107.86820, 1.44 ,  8]
    ptable[ 48] = ['Cd', 112.41400, 1.49 ,  8]
    ptable[ 49] = ['In', 114.81800, 1.66 ,  8]
    ptable[ 50] = ['Sn', 118.71000, 1.62 ,  8]
    ptable[ 51] = ['Sb', 121.76000, 1.41 ,  8]
    ptable[ 52] = ['Te', 127.60000, 1.37 ,  8]
    ptable[ 53] = [ 'I', 126.90447, 1.33 ,  1]
    ptable[ 54] = ['Xe', 131.29300, 2.17 ,  0]
    ptable[ 55] = ['Cs', 132.90545, 2.62 ,  8]
    ptable[ 56] = ['Ba', 137.32700, 2.17 ,  8]
    ptable[ 57] = ['La', 138.90547, 1.88 ,  8]
    ptable[ 58] = ['Ce', 140.11600, 1.818,  8]
    ptable[ 59] = ['Pr', 140.90766, 1.824,  8]
    ptable[ 60] = ['Nd', 144.24200, 1.814,  8]
    ptable[ 61] = ['Pm', 145.00000, 1.834,  8]
    ptable[ 62] = ['Sm', 150.36000, 1.804,  8]
    ptable[ 63] = ['Eu', 151.96400, 2.084,  8]
    ptable[ 64] = ['Gd', 157.25000, 1.804,  8]
    ptable[ 65] = ['Tb', 158.92535, 1.773,  8]
    ptable[ 66] = ['Dy', 162.50000, 1.781,  8]
    ptable[ 67] = ['Ho', 164.93033, 1.762,  8]
    ptable[ 68] = ['Er', 167.25900, 1.761,  8]
    ptable[ 69] = ['Tm', 168.93422, 1.759,  8]
    ptable[ 70] = ['Yb', 173.04500, 1.922,  8]
    ptable[ 71] = ['Lu', 174.96680, 1.738,  8]
    ptable[ 72] = ['Hf', 178.49000, 1.57 ,  8]
    ptable[ 73] = ['Ta', 180.94788, 1.43 ,  8]
    ptable[ 74] = [ 'W', 183.84000, 1.37 ,  8]
    ptable[ 75] = ['Re', 186.20700, 1.37 ,  8]
    ptable[ 76] = ['Os', 190.23000, 1.34 ,  8]
    ptable[ 77] = ['Ir', 192.21700, 1.35 ,  8]
    ptable[ 78] = ['Pt', 195.08400, 1.38 ,  8]
    ptable[ 79] = ['Au', 196.96657, 1.44 ,  8]
    ptable[ 80] = ['Hg', 200.59200, 1.52 ,  8]
    ptable[ 81] = ['Tl', 204.38350, 1.71 ,  8]
    ptable[ 82] = ['Pb', 207.20000, 1.75 ,  8]
    ptable[ 83] = ['Bi', 208.98040, 1.70 ,  8]
    ptable[ 84] = ['Po', 209.00000, 1.40 ,  8]
    ptable[ 85] = ['At', 210.00000, 1.40 ,  1]
    ptable[ 86] = ['Rn', 222.00000, 2.40 ,  0]
    ptable[ 87] = ['Fr', 223.00000, 2.70 ,  8]
    ptable[ 88] = ['Ra', 226.00000, 2.20 ,  8]
    ptable[ 89] = ['Ac', 227.00000, 2.00 ,  8]
    ptable[ 90] = ['Th', 232.03770, 1.79 ,  8]
    ptable[ 91] = ['Pa', 231.03588, 1.63 ,  8]
    ptable[ 92] = [ 'U', 238.02891, 1.56 ,  8]
    ptable[ 93] = ['Np', 237.00000, 1.55 ,  8]
    ptable[ 94] = ['Pu', 244.00000, 1.59 ,  8]
    ptable[ 95] = ['Am', 243.00000, 1.73 ,  8]
    ptable[ 96] = ['Cm', 247.00000, 1.74 ,  8]
    ptable[ 97] = ['Bk', 247.00000, 1.70 ,  8]
    ptable[ 98] = ['Cf', 251.00000, 1.86 ,  8]
    ptable[ 99] = ['Es', 252.00000, 1.86 ,  8]
    ptable[100] = ['Fm', 257.00000, 2.00 ,  8]
    ptable[101] = ['Md', 258.00000, 2.00 ,  8]
    ptable[102] = ['No', 259.00000, 2.00 ,  8]
    ptable[103] = ['Lr', 266.00000, 2.00 ,  8]
    ptable[104] = ['Rf', 267.00000, 2.00 ,  8]
    ptable[105] = ['Db', 268.00000, 2.00 ,  8]
    ptable[106] = ['Sg', 269.00000, 2.00 ,  8]
    ptable[107] = ['Bh', 270.00000, 2.00 ,  8]
    ptable[108] = ['Hs', 277.00000, 2.00 ,  8]
    ptable[109] = ['Mt', 278.00000, 2.00 ,  8]
    ptable[110] = ['Ds', 281.00000, 2.00 ,  8]
    ptable[111] = ['Rg', 282.00000, 2.00 ,  8]
    ptable[112] = ['Cn', 285.00000, 2.00 ,  8]
    ptable[113] = ['Nh', 286.00000, 2.00 ,  8]
    ptable[114] = ['Fl', 289.00000, 2.00 ,  8]
    ptable[115] = ['Mc', 290.00000, 2.00 ,  8]
    ptable[116] = ['Lv', 293.00000, 2.00 ,  8]
    ptable[117] = ['Ts', 294.00000, 2.00 ,  8]
    ptable[118] = ['Og', 294.00000, 2.00 ,  8]

    def __init__(self):
        raise PTError('The PeriodicTable class cannot be instantiated.')

    @classmethod
    def get_atnum(cls, symbol):
        """Converts symbol to atomic number"""
        for key in cls.ptable.keys():
            if cls.ptable[key][0] == symbol.capitalize():
                return key

    @classmethod
    def get_symbol(cls, atnum):
        """Converts atomic number to symbol"""
        return cls.ptable[atnum][0]
    
    @classmethod
    def get_radius(cls, symbol):
        """Returns atomic radius for a given element"""
        atnum = cls.get_atnum(symbol)
        return float(cls.ptable[atnum][2])
    
    @classmethod
    def get_mass(cls, symbol):
        """Returns atomic mass for a given element"""
        atnum = cls.get_atnum(symbol)
        return float(cls.ptable[atnum][1])

    @classmethod
    def get_connectors(cls, symbol):
        """Returns number of possible attachments to a given element"""
        atnum = cls.get_atnum(symbol)
        return int(cls.ptable[atnum][-1])
