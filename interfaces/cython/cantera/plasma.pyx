# This file is part of Cantera. See License.txt in the top-level directory or
# at https://www.cantera.org/license.txt for license and copyright information.

import warnings
import weakref

cdef class PlasmaPhase(ThermoPhase):
    """ A class representing a plasma phase"""
    def __cinit__(self, *args, **kwargs):
        if pystr(self.thermo.type()) != "Plasma":
            raise TypeError('Underlying ThermoPhase object is of the wrong type.')
        self.plasma = <CxxPlasmaPhase*>(self.thermo)

    property electron_energy_grid:
        """ Electron energy grid [eV]"""
        def __get__(self):
            cdef vector[double] cxxdata
            self.plasma.getElectronEnergyGrid(cxxdata)
            return np.fromiter(cxxdata, np.double)
        def __set__(self, grid):
            cdef vector[double] cxxdata
            for value in grid:
                cxxdata.push_back(value)
            self.plasma.setElectronEnergyGrid(cxxdata)
