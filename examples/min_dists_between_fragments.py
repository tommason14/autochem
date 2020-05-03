from autochem import Molecule, df_from_namedtuples
from glob import glob
import itertools
from collections import namedtuple

info = namedtuple("info", "xyz atom1 mol1 atom2 mol2 dist")

def distances(namedtup, molecule):
    return [
        namedtup(
            molecule.xyz,
            f"{atom1.symbol}_{atom1.index}",
            atom1.fragment,
            f"{atom2.symbol}_{atom2.index}",
            atom2.fragment,
            atom1.distance_to(atom2),
        )
        for frag1, frag2 in itertools.combinations(molecule.fragments.values(), 2)
        for atom1 in frag1["atoms"]
        for atom2 in frag2["atoms"]
    ]

dists = []
for xyz in glob('*xyz'):
    mol = Molecule(xyz)
    dists += distances(info, mol)

df = df_from_namedtuples(info, dists)
mindists = df.groupby(["xyz", "mol1", "mol2"])["dist"].min().reset_index().sort_values(
    ["xyz", "mol1", "mol2"]
)
mindists.to_csv('min_dists_between_frags.csv', index=False)
