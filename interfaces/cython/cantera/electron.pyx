# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

import warnings
import weakref

cdef np.ndarray get_electron_1d(Electron elect, electronMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(elect.thermo.nSpecies())
    method(elect.electron, &data[0])
    if elect._selected_species.size:
        return data[elect._selected_species]
    else:
        return data

cdef class Electron(_SolutionBase):
    """
    This class is used to compute electron properties for a phase of matter.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def set_electron_energy_grid(self, grid):
        """ Set the grid of cell boundary of electron energy [eV]"""
        cdef np.ndarray[np.double_t, ndim=1] data = \
            np.ascontiguousarray(grid, dtype=np.double)
        self.electron.setupGrid(len(data), &data[0])

    property electron_temperature:
        """electron temperature"""
        def __get__(self):
            return self.electron.electronTemperature()

    property electron_mobility:
        """electron mobility [m^2/V/s)]"""
        def __get__(self):
            return self.electron.electronMobility()

    property electron_diffusivity:
        """electron diffusivity [m^2/s]"""
        def __get__(self):
            return self.electron.electronDiffusivity()

    def electron_rate_coefficient(self, k):
        """rate coefficient of process k"""
        return self.electron.rateCoefficient(k)

    def electron_inverse_rate_coefficient(self, k):
        """inverse rate coefficient of process k"""
        return self.electron.inverseRateCoefficient(k)

    property electron_power_gain:
        """
        Electron power gain. [eV/s]
        """
        def __get__(self):
            return self.electron.powerGain()

    property electron_elastic_power_loss:
        """
        Electron elastic power loss. [eV/s]
        """
        def __get__(self):
            return self.electron.elasticPowerLoss()

    property electron_inelastic_power_loss:
        """
        Electron inelastic power loss. [eV/s]
        """
        def __get__(self):
            return self.electron.inelasticPowerLoss()

    property net_plasma_production_rates:
        """
        Net plasma production rates for each species. [kmol/m^3/s] for bulk phases or
        [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_electron_1d(self, elect_getNetPlasmaProductionRates)

    def set_chemionization_scattering_rate(self, rate):
        """ Set chemionization scattering-in rate """
        self.electron.setChemionScatRate(rate)

    def set_boltzmann_solver(self, maxn=100, rtol=1e-5, delta0=1e14,
                             m=4.0, init_kTe=0.0, warn=True):
        """ Set boltzmann solver"""
        self.electron.setBoltzmannSolver(maxn, rtol, delta0, m, init_kTe, <cbool>warn)

    property mean_electron_energy:
        """mean electron energy [eV]"""
        def __get__(self):
            return self.electron.meanElectronEnergy()

    property electric_field:
        """reduced electric field strength [V-m2]"""
        def __get__(self):
            return self.electron.electricField()
        def __set__(self, E):
            self.electron.setElectricField(E)

    property electric_field_freq:
        """electric field freq [Hz]"""
        def __get__(self):
            return self.electron.electricFieldFreq()
        def __set__(self, F):
            self.electron.setElectricFieldFreq(F)


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
