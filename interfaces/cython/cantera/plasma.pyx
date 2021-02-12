# This file is part of Cantera. See License.txt in the top-level directory or
# at https://www.cantera.org/license.txt for license and copyright information.

import warnings
import weakref

cdef class PlasmaPhase(ThermoPhase):
    """ A class representing a plasma phase"""
    def __cinit__(self, *args, **kwargs):
        if pystr(self.thermo.type()) != "WeaklyIonizedGas":
            raise TypeError('Underlying ThermoPhase object is of the wrong type.')
        self.plasma = <CxxPlasmaPhase*>(self.thermo)

    def set_electron_energy_grid(self, grid):
        """ Set the grid of cell boundary of electron energy [eV]"""
        cdef np.ndarray[np.double_t, ndim=1] data = \
            np.ascontiguousarray(grid, dtype=np.double)
        self.plasma.setupGrid(len(data), &data[0])

    property electron_temperature:
        """electron temperature"""
        def __get__(self):
            return self.plasma.electronTemperature()

    property electron_mobility:
        """electron mobility [m^2/V/s]"""
        def __get__(self):
            return self.plasma.electronMobility()

    property electron_diffusivity:
        """electron diffusivity [m^2/s]"""
        def __get__(self):
            return self.plasma.electronDiffusivity()

    property electron_total_collision_frequency:
        """electron total collision frequency"""
        def __get__(self):
            return self.plasma.totalCollisionFreq()

    def plasma_process_rate_coefficient(self, k):
        """
        rate coefficient of process k in inverse
        number density [m^3/s].
        """
        return self.plasma.rateCoefficient(k)

    def plasma_process_reverse_rate_coefficient(self, k):
        """
        reverse rate coefficient of process k in inverse
        number density [m^3/s].
        """
        return self.plasma.reverseRateCoefficient(k)

    property electron_power_gain:
        """
        Electron power gain. [eV/s]
        """
        def __get__(self):
            return self.plasma.powerGain()

    property electron_elastic_power_loss:
        """
        The elastic power loss of one electron. [eV/s]
        """
        def __get__(self):
            return self.plasma.elasticPowerLoss()

    property electron_inelastic_power_loss:
        """
        The inelastic power loss of one electron. [eV/s]
        """
        def __get__(self):
            return self.plasma.inelasticPowerLoss()

    property electron_total_power_loss:
        """
        The total power loss of one electron. [eV/s]
        """
        def __get__(self):
            return (self.electron_elastic_power_loss +
                    self.electron_inelastic_power_loss)

    def setup_boltzmann_solver(self, maxn=100, rtol=1e-5,
                               delta0=1e14, m=4.0):
        """ 
        Setup Boltzmann Equation solver.
        :param maxn:
            Maximum number of iterations
        :param rtol:
            Relative tolerance
        :param delta0:
            Initial value of the iteration parameter
        :param m:
            Reduction factor of error
        """
        self.plasma.setupBoltzmannSolver(maxn, rtol, delta0, m)

    def set_mole_fraction_threshold(self, fraction):
        """ Set mole fraction threshold. """
        self.plasma.setMoleFractionThreshold(fraction)

    def set_initial_mean_electron_energy(self, init_kTe):
        """ Set initial mean electron energy. [eV]"""
        self.plasma.setInitialMeanElectronEnergy(init_kTe)

    def set_reuse_EEDF(self, reuse):
        """ Set flag of reusing old EEDF as initial EEDF"""
        self.plasma.setReuseEEDF(<cbool>reuse)

    property mean_electron_energy:
        """mean electron energy [eV]"""
        def __get__(self):
            return self.plasma.meanElectronEnergy()

    property electric_field:
        """reduced electric field strength [V-m2]"""
        def __get__(self):
            return self.plasma.electricField()
        def __set__(self, E):
            self.plasma.setElectricField(E)

    property electric_field_freq:
        """electric field freq [Hz]"""
        def __get__(self):
            return self.plasma.electricFieldFreq()
        def __set__(self, F):
            self.plasma.setElectricFieldFreq(F)

    def electron_collision_target(self, k):
        """ The target of the collision process. """
        return pystr(self.plasma.target(k))

    def electron_collision_kind(self, k):
        """ The kind of the collision process. """
        return pystr(self.plasma.kind(k))

    def electron_collision_product(self,k):
        """ The product of the collision. process """
        return pystr(self.plasma.product(k))

    def electron_collision_threshold(self, k):
        """ The threhols of the collision process. """
        return self.plasma.threshold(k)


cdef class ElectronCrossSection:
    """
    A class which stores data about a cross section of electron.

    :param kind:
        A string giving the kind of the collision, e.g. ``'EFFECTIVE'``
    :param target:
        A string giving the target of the collision, e.g. ``'N2'``
    :param init:
        Used internally when wrapping :ct:`ElectronCrossSection` objects returned from C++

    Example:
    S = ct.ElectronCrossSection.listFromFile('gri30.yaml')
    """
    def __cinit__(self, *args, init=True, **kwargs):
        self._electron_cross_section.reset(new CxxElectronCrossSection())
        self.electron_cross_section = self._electron_cross_section.get()

    def __init__(self, kind=None, target=None, data=None,
                 *args, init=True, **kwargs):
        if not init:
            return

        if kind is not None:
            self.electron_cross_section.kind = stringify(kind)

        if target is not None:
            self.electron_cross_section.target = stringify(target)

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
        if filename.endswith('.yml') or filename.endswith('.yaml'):
            root = AnyMapFromYamlFile(stringify(filename))
            cxx_electron_cross_section = CxxGetElectronCrossSection(root[stringify(section)])
        else:
            raise ValueError("The file type must be yml or yaml")

        cross_section = []
        for a in cxx_electron_cross_section:
            b = ElectronCrossSection(init=False)
            b._assign(a)
            cross_section.append(b)
        return cross_section

    property kind:
        """ The kind of the collision process. """
        def __get__(self):
            return pystr(self.electron_cross_section.kind)

    property target:
        """ The target of the collision process. """
        def __get__(self):
            return pystr(self.electron_cross_section.target)

    property product:
        """ The product of the collision process. """
        def __get__(self):
            return pystr(self.electron_cross_section.product)

    property threshold:
        """ The threhols of the collision process. """
        def __get__(self):
            return self.electron_cross_section.threshold

    property process:
        """ The process of the collision. 
            Ex. "O2 => O2^+" (ionization)
                "O2 => O + O" (dissociation)
                "O2 => O2(a1)" (exitation)"""
        def __get__(self):
            return self.target + " => " + self.product

    property data:
        """ The data of electron collision cross section. """
        def __get__(self):
            return self.electron_cross_section.data
