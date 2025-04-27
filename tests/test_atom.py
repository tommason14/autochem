import autochem
import re


def test_atom():
    atom = autochem.Atom(symbol="H", coords=[1, 2, 3])
    assert re.sub(r"\s+", " ", str(atom)) == "Atom: H 1.00000 2.00000 3.00000"
