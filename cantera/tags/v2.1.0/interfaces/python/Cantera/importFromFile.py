"""Functions to import phase and interface definitions from CTI or
CTML files. This module is imported when the Cantera package is
imported, and therfore does not need to be explicitly imported in
application programs."""

import solution
import Interface
import Edge
import XML

def importPhase(file, name = '', loglevel = 0, debug = 0):
    """Import one phase from an input  file.  If 'name' is specified, the
    phase definition with this name will be imported, otherwise the first
    phase definition in the file will be imported. If 'loglevel' is set to
    a positive integer, additional information will be printed or written
    to log files about the details of the object constructed.
    """
    return importPhases(file, [name], loglevel, debug)[0]

def importPhases(file, names = [], loglevel = 0, debug = 0):
    """Import multiple phases from one file. The phase names should be
    entered as a list of strings. See: importPhase """
    s = []
    for nm in names:
        s.append(solution.Solution(src=file,id=nm,loglevel=loglevel,debug=debug))
    return s

def importInterface(file, name = '', phases = []):
    """Import an interface definition from input file 'file', and return
    an instance of class Interface implementing this definition. If 'name'
    is specified, the definition by this name will be imported; otherwise, the
    first interface definition in the file will be imported.

    The 'phases' argument is a list of objects representing the other phases
    that participate in the interfacial reactions, for example an object
    representing a gas phase or a solid.

    >>> gas1, cryst1 = importPhases('diamond.cti', ['gas', 'solid'])
    >>> diamond_surf = importInterface('diamond.cti', [gas1, cryst1])

    Note the difference between the lists in the argument lists of these
    two functions. In importPhases, a list of name strings is entered,
    which are used to identify the appropriate definitions in the input
    file to build the objects. In importInterface, the list is of the
    objects that were built by importPhases. The reason these objects must be
    given as inputs is that these objects will be queried when phase
    properties (temperature, pressure, composition,
    electric potential) are needed to compute the reaction rates of progress.
    """
    if name:
        src = file+'#'+name
    else:
        src = file
    return Interface.Interface(src = src, phases = phases)


def importEdge(file, name = '', surfaces = []):
    if name:
        src = file+'#'+name
    else:
        src = file
    return Edge.Edge(src = src, surfaces = surfaces)
