from . import Solution


def import_phases(filename, phase_names):
    """
    Import multiple phases from one file. The phase names should be
    entered as a list of strings.
    """
    return [Solution(filename, p) for p in phase_names]
