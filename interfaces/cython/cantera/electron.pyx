# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

import warnings
import weakref


cdef class ElectronCrossSection:
    """
    A class which stores data about a cross section of electron.

    :param kind:
        A string giving the kind of the collision, e.g. ``'EFFECTIVE'``
    :param target:
        A string giving the target of the collision, e.g. ``'N2'``
    :param mass_ratio:
        The mass ratio of electron to molecule.
    :param init:
        Used internally when wrapping :ct:`ElectronCrossSection` objects returned from C++

    Example:
    S = ct.ElectronCrossSection.listFromFile('gri30.yaml')
    """
    def __cinit__(self, *args, init=True, **kwargs):
        self._electron_cross_section.reset(new CxxElectronCrossSection())
        self.electron_cross_section = self._electron_cross_section.get()

    def __init__(self, kind=None, target=None, mass_ratio=None, data=None,
                 *args, init=True, **kwargs):
        if not init:
            return

        if kind is not None:
            self.electron_cross_section.kind = stringify(kind)

        if target is not None:
            self.electron_cross_section.target = stringify(target)

        if mass_ratio is not None:
            self.electron_cross_section.mass_ratio = mass_ratio

        if data is not None:
            self.electron_cross_section.data = data

    cdef _assign(self, shared_ptr[CxxElectronCrossSection] other):
        self._electron_cross_section = other
        self.electron_cross_section = self._electron_cross_section.get()

    @staticmethod
    def listFromFile(filename, section='cross_section'):
        """
        Create a list of ElectronCrossSection objects from all of the
        cross section defined in a YAML.

        Directories on Cantera's input file path will be searched for the
        specified file.
        """
        if filename.lower().split('.')[-1] in ('yml', 'yaml'):
            root = AnyMapFromYamlFile(stringify(filename))
            cxx_electron_cross_section = CxxGetElectronCrossSection(root[stringify(section)])

        cross_section = []
        for a in cxx_electron_cross_section:
            b = ElectronCrossSection(init=False)
            b._assign(a)
            cross_section.append(b)
        return cross_section

    property kind:
        """ The kind of the collision. """
        def __get__(self):
            return pystr(self.electron_cross_section.kind)

    property target:
        """ The target of the collision. """
        def __get__(self):
            return pystr(self.electron_cross_section.target)

    property mass_ratio:
        """ The mass ratio of electron to molecule. """
        def __get__(self):
            return self.electron_cross_section.mass_ratio

    property data:
        """ The mass ratio of electron to molecule. """
        def __get__(self):
            return self.electron_cross_section.data
