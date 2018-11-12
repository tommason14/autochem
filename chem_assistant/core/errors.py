__all__ = ['PlamsError', 'FileError', 'ResultsError', 'JobError', 'PTError', 'MoleculeError', 'SuperCompError']

class PlamsError(Exception):
    """General PLAMS error."""
    pass

class FileError(Exception):
    """File or filesystem related error."""
    pass

class ResultsError(Exception):
    """|Results| related error."""
    pass

class JobError(Exception):
    """|Job| related error."""
    pass

class PTError(Exception):
    """|PeriodicTable| error."""
    pass

class MoleculeError(Exception):
    """|Molecule| related error."""
    pass

class SuperCompError(Exception):
    pass