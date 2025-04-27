from autochem.core.periodic_table import PeriodicTable
import math


def test_periodic_table():
    assert PeriodicTable.get_atomic_number("H") == 1
    assert math.ceil(PeriodicTable.get_mass("O")) == 16
