__all__ = ["PeriodicTable"]

from collections import namedtuple

# no element included here, use that as the key to the dict
atom = namedtuple("atom", "atomic_number mass radius connectors vdw_radius")


class PeriodicTable:
    """Helper class to allow for lookup of atomic properties. Can convert between symbol and atomic number"""

    ptable = {}
    # symbol =  atomic number, radius, connectors, vdw radii
    # atomic weights from: http://www.ciaaw.org/atomic-weights.htm
    # need to add more vdw radii from molden oglmol file (vdwr array)

    ptable["Xx"] = atom(0, 0.00000, 0.00, 0, 0.000)
    ptable["H"] = atom(1, 1.00798, 0.30, 1, 0.430)
    ptable["He"] = atom(2, 4.00260, 0.99, 0, 0.741)
    ptable["Li"] = atom(3, 6.96750, 1.52, 8, 0.880)
    ptable["Be"] = atom(4, 9.01218, 1.12, 8, 0.550)
    ptable["B"] = atom(5, 10.81350, 0.88, 6, 1.030)
    ptable["C"] = atom(6, 12.01060, 0.77, 4, 0.900)
    ptable["N"] = atom(7, 14.00685, 0.70, 3, 0.880)
    ptable["O"] = atom(8, 15.99940, 0.66, 2, 0.880)
    ptable["F"] = atom(9, 18.99840, 0.64, 1, 0.840)
    ptable["Ne"] = atom(10, 20.17970, 1.60, 0, 0.815)
    ptable["Na"] = atom(11, 22.98977, 1.86, 8, 1.170)
    ptable["Mg"] = atom(12, 24.30550, 1.60, 8, 1.300)
    ptable["Al"] = atom(13, 26.98154, 1.43, 8, 1.550)
    ptable["Si"] = atom(14, 28.08500, 1.17, 8, 1.400)
    ptable["P"] = atom(15, 30.97376, 1.10, 8, 1.250)
    ptable["S"] = atom(16, 32.06750, 1.04, 2, 1.220)
    ptable["Cl"] = atom(17, 35.45150, 0.99, 1, 1.190)
    ptable["Ar"] = atom(18, 39.94800, 1.92, 0, 0.995)
    ptable["K"] = atom(19, 39.09830, 2.31, 8, 1.530)
    ptable["Ca"] = atom(20, 40.07800, 1.97, 8, 1.190)
    ptable["Sc"] = atom(21, 44.95591, 1.60, 8, 1.640)
    ptable["Ti"] = atom(22, 47.86700, 1.46, 8, 1.670)
    ptable["V"] = atom(23, 50.94150, 1.31, 8, 1.530)
    ptable["Cr"] = atom(24, 51.99610, 1.25, 8, 1.550)
    ptable["Mn"] = atom(25, 54.93804, 1.29, 8, 1.555)
    ptable["Fe"] = atom(26, 55.84500, 1.26, 8, 1.540)
    ptable["Co"] = atom(27, 58.93319, 1.25, 8, 1.530)
    ptable["Ni"] = atom(28, 58.69340, 1.24, 8, 1.700)
    ptable["Cu"] = atom(29, 63.54600, 1.28, 8, 1.720)
    ptable["Zn"] = atom(30, 65.38000, 1.33, 8, 1.650)
    ptable["Ga"] = atom(31, 69.72300, 1.41, 8, 1.420)
    ptable["Ge"] = atom(32, 72.63000, 1.22, 8, 1.370)
    ptable["As"] = atom(33, 74.92159, 1.21, 8, 1.410)
    ptable["Se"] = atom(34, 78.97100, 1.17, 8, 1.420)
    ptable["Br"] = atom(35, 79.90400, 1.14, 1, 1.410)
    ptable["Kr"] = atom(36, 83.79800, 1.97, 0, 1.069)
    ptable["Rb"] = atom(37, 85.46780, 2.44, 8, 1.670)
    ptable["Sr"] = atom(38, 87.62000, 2.15, 8, 1.320)
    ptable["Y"] = atom(39, 88.90584, 1.80, 8, 1.980)
    ptable["Zr"] = atom(40, 91.22400, 1.57, 8, 1.760)
    ptable["Nb"] = atom(41, 92.90637, 1.41, 8, 1.680)
    ptable["Mo"] = atom(42, 95.95000, 1.36, 8, 1.670)
    ptable["Tc"] = atom(43, 98.00000, 1.35, 8, 1.550)
    ptable["Ru"] = atom(44, 101.07000, 1.33, 8, 1.600)
    ptable["Rh"] = atom(45, 102.90550, 1.34, 8, 1.650)
    ptable["Pd"] = atom(46, 106.42000, 1.38, 8, 1.700)
    ptable["Ag"] = atom(47, 107.86820, 1.44, 8, 1.790)
    ptable["Cd"] = atom(48, 112.41400, 1.49, 8, 1.890)
    ptable["In"] = atom(49, 114.81800, 1.66, 8, 1.830)
    ptable["Sn"] = atom(50, 118.71000, 1.62, 8, 1.660)
    ptable["Sb"] = atom(51, 121.76000, 1.41, 8, 0.000)
    ptable["Te"] = atom(52, 127.60000, 1.37, 8, 0.000)
    ptable["I"] = atom(53, 126.90447, 1.33, 1, 0.000)
    ptable["Xe"] = atom(54, 131.29300, 2.17, 0, 0.000)
    ptable["Cs"] = atom(55, 132.90545, 2.62, 8, 0.000)
    ptable["Ba"] = atom(56, 137.32700, 2.17, 8, 0.000)
    ptable["La"] = atom(57, 138.90547, 1.88, 8, 0.000)
    ptable["Ce"] = atom(58, 140.11600, 1.818, 8, 0.000)
    ptable["Pr"] = atom(59, 140.90766, 1.824, 8, 0.000)
    ptable["Nd"] = atom(60, 144.24200, 1.814, 8, 0.000)
    ptable["Pm"] = atom(61, 145.00000, 1.834, 8, 0.000)
    ptable["Sm"] = atom(62, 150.36000, 1.804, 8, 0.000)
    ptable["Eu"] = atom(63, 151.96400, 2.084, 8, 0.000)
    ptable["Gd"] = atom(64, 157.25000, 1.804, 8, 0.000)
    ptable["Tb"] = atom(65, 158.92535, 1.773, 8, 0.000)
    ptable["Dy"] = atom(66, 162.50000, 1.781, 8, 0.000)
    ptable["Ho"] = atom(67, 164.93033, 1.762, 8, 0.000)
    ptable["Er"] = atom(68, 167.25900, 1.761, 8, 0.000)
    ptable["Tm"] = atom(69, 168.93422, 1.759, 8, 0.000)
    ptable["Yb"] = atom(70, 173.04500, 1.922, 8, 0.000)
    ptable["Lu"] = atom(71, 174.96680, 1.738, 8, 0.000)
    ptable["Hf"] = atom(72, 178.49000, 1.57, 8, 0.000)
    ptable["Ta"] = atom(73, 180.94788, 1.43, 8, 0.000)
    ptable["W"] = atom(74, 183.84000, 1.37, 8, 0.000)
    ptable["Re"] = atom(75, 186.20700, 1.37, 8, 0.000)
    ptable["Os"] = atom(76, 190.23000, 1.34, 8, 0.000)
    ptable["Ir"] = atom(77, 192.21700, 1.35, 8, 0.000)
    ptable["Pt"] = atom(78, 195.08400, 1.38, 8, 0.000)
    ptable["Au"] = atom(79, 196.96657, 1.44, 8, 0.000)
    ptable["Hg"] = atom(80, 200.59200, 1.52, 8, 0.000)
    ptable["Tl"] = atom(81, 204.38350, 1.71, 8, 0.000)
    ptable["Pb"] = atom(82, 207.20000, 1.75, 8, 0.000)
    ptable["Bi"] = atom(83, 208.98040, 1.70, 8, 0.000)
    ptable["Po"] = atom(84, 209.00000, 1.40, 8, 0.000)
    ptable["At"] = atom(85, 210.00000, 1.40, 1, 0.000)
    ptable["Rn"] = atom(86, 222.00000, 2.40, 0, 0.000)
    ptable["Fr"] = atom(87, 223.00000, 2.70, 8, 0.000)
    ptable["Ra"] = atom(88, 226.00000, 2.20, 8, 0.000)
    ptable["Ac"] = atom(89, 227.00000, 2.00, 8, 0.000)
    ptable["Th"] = atom(90, 232.03770, 1.79, 8, 0.000)
    ptable["Pa"] = atom(91, 231.03588, 1.63, 8, 0.000)
    ptable["U"] = atom(92, 238.02891, 1.56, 8, 0.000)
    ptable["Np"] = atom(93, 237.00000, 1.55, 8, 0.000)
    ptable["Pu"] = atom(94, 244.00000, 1.59, 8, 0.000)
    ptable["Am"] = atom(95, 243.00000, 1.73, 8, 0.000)
    ptable["Cm"] = atom(96, 247.00000, 1.74, 8, 0.000)
    ptable["Bk"] = atom(97, 247.00000, 1.70, 8, 0.000)
    ptable["Cf"] = atom(98, 251.00000, 1.86, 8, 0.000)
    ptable["Es"] = atom(99, 252.00000, 1.86, 8, 0.000)
    ptable["Fm"] = atom(100, 257.00000, 2.00, 8, 0.000)
    ptable["Md"] = atom(101, 258.00000, 2.00, 8, 0.000)
    ptable["No"] = atom(102, 259.00000, 2.00, 8, 0.000)
    ptable["Lr"] = atom(103, 266.00000, 2.00, 8, 0.000)
    ptable["Rf"] = atom(104, 267.00000, 2.00, 8, 0.000)
    ptable["Db"] = atom(105, 268.00000, 2.00, 8, 0.000)
    ptable["Sg"] = atom(106, 269.00000, 2.00, 8, 0.000)
    ptable["Bh"] = atom(107, 270.00000, 2.00, 8, 0.000)
    ptable["Hs"] = atom(108, 277.00000, 2.00, 8, 0.000)
    ptable["Mt"] = atom(109, 278.00000, 2.00, 8, 0.000)
    ptable["Ds"] = atom(110, 281.00000, 2.00, 8, 0.000)
    ptable["Rg"] = atom(111, 282.00000, 2.00, 8, 0.000)
    ptable["Cn"] = atom(112, 285.00000, 2.00, 8, 0.000)
    ptable["Nh"] = atom(113, 286.00000, 2.00, 8, 0.000)
    ptable["Fl"] = atom(114, 289.00000, 2.00, 8, 0.000)
    ptable["Mc"] = atom(115, 290.00000, 2.00, 8, 0.000)
    ptable["Lv"] = atom(116, 293.00000, 2.00, 8, 0.000)
    ptable["Ts"] = atom(117, 294.00000, 2.00, 8, 0.000)
    ptable["Og"] = atom(118, 294.00000, 2.00, 8, 0.000)

    def __init__(self):
        raise AttributeError("The PeriodicTable class cannot be instantiated.")

    @classmethod
    def get_atomic_number(cls, symbol):
        """Converts symbol to atomic number"""
        return cls.ptable.get(symbol).atomic_number

    @classmethod
    def get_symbol(cls, atomic_number):
        """Converts atomic number to symbol"""
        for symbol, atoms in cls.ptable.items():
            if atom.atomic_number == atomic_number:
                return symbol

    @classmethod
    def get_radius(cls, symbol):
        """Returns atomic radius for a given element"""
        return cls.ptable.get(symbol).radius

    @classmethod
    def get_mass(cls, symbol):
        """Returns atomic mass for a given element"""
        return cls.ptable.get(symbol).mass

    @classmethod
    def get_connectors(cls, symbol):
        """Returns number of possible attachments to a given element"""
        return cls.ptable.get(symbol).connectors

    @classmethod
    def get_vdw(cls, symbol):
        """Returns van der waals radius of a given element"""
        return cls.ptable.get(symbol).vdw_radius
