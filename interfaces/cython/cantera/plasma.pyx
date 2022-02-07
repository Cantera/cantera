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

    def set_electron_energy_distribution(self, grid, distrb):
        """ Set electron energy distribution. When this method is used, electron
        temeprature is calculated from the distribution.
        :param grid:
            vector of electron energy grid [eV]
        :param distrb:
            vector of distribution
        """
        cdef vector[double] cxxdata_grid
        cdef vector[double] cxxdata_distrb
        for value in grid:
            cxxdata_grid.push_back(value)
        for value in distrb:
            cxxdata_distrb.push_back(value)
        self.plasma.setElectronEnergyDistrb(cxxdata_grid, cxxdata_distrb)

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

    property electron_energy_distribution:
        """ Electron energy distribution """
        def __get__(self):
            cdef vector[double] cxxdata
            self.plasma.getElectronEnergyDistrb(cxxdata)
            return np.fromiter(cxxdata, np.double)

    property mean_electron_energy:
        """ Mean electron energy [eV] """
        def __get__(self):
            return self.plasma.meanElectronEnergy()
