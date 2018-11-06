__all__ = ['PlamsError', 'FileError', 'ResultsError', 'JobError', 'PTError', 'UnitsError', 'MoleculeError']

class PlamsError(Exception):
    """General PLAMS error."""
    pass

class FileError(PlamsError):
    """File or filesystem related error."""
    pass

class ResultsError(PlamsError):
    """|Results| related error."""
    pass

class JobError(PlamsError):
    """|Job| related error."""
    pass

class PTError(PlamsError):
    """:class:`Periodic table<scm.plams.utils.PeriodicTable>` error."""
    pass

class UnitsError(PlamsError):
    """:class:`Units converter<scm.plams.utils.Units>` error."""
    pass

class MoleculeError(PlamsError):
    """|Molecule| related error."""
    pass
