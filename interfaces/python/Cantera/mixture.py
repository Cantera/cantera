"""
Multiphase mixtures.
"""

import _cantera
import types
from Cantera.num import zeros, array, asarray
from exceptions import CanteraError
from Cantera import writeLogFile

class Mixture:
    """
    Multiphase mixtures.  Class Mixture represents
    mixtures of one or more phases of matter.  To construct a mixture,
    supply a list of phases to the constructor, each paired with the
    number of moles for that phase:

    >>> gas = importPhase('gas.cti')
    >>> gas.speciesNames()
    ['H2', 'H', 'O2', 'O', 'OH']
    >>> graphite = importPhase('graphite.cti')
    >>> graphite.speciesNames()
    ['C(g)']
    >>> mix = Mixture([(gas, 1.0), (graphite, 0.1)])
    >>> mix.speciesNames()
    ['H2', 'H', 'O2', 'O', 'OH', 'C(g)']

    Note that the objects representing each phase compute only the
    intensive state of the phase -- they do not store any information
    on the amount of this phase. Mixture objects, on the other hand, represent
    the full extensive state.

    Mixture objects are 'lightweight' in the sense that they do not
    store parameters needed to compute thermodynamic or kinetic
    properties of the phases. These are contained in the
    ('heavyweight') phase objects. Multiple mixture objects may be
    constructed using the same set of phase objects. Each one stores
    its own state information locally, and synchronizes the phases
    objects whenever it requires phase properties.

    """

    def __init__(self, phases=[]):
        self.__mixid = _cantera.mix_new()
        self._spnames = []
        self._phases = []
        if phases:
            for p in phases:
                try:
                    ph = p[0]
                    moles = p[1]
                except:
                    ph = p
                    if p == phases[0]:
                        moles = 1
                    else:
                        moles = 0
                self._addPhase(ph, moles)
                self._phases.append(ph)
        _cantera.mix_init(self.__mixid)
        self.setTemperature(self._phases[0].temperature())
        self.setPressure(self._phases[0].pressure())

    def __del__(self):
        """Delete the Mixture instance. The phase objects are not deleted."""
        _cantera.mix_del(self.__mixid)

    def __str__(self):
        s = ''
        for p in range(len(self._phases)):
            s += '\n*******************    Phase '+self._phases[p].name()+'    ******************************\n'
            s += '\n Moles: '+`self.phaseMoles(p)`+'\n'
            s += self._phases[p].__repr__()+'\n\n'
        return s

    def _addPhase(self, phase = None, moles = 0.0):
        """Add a phase to the mixture."""
        for k in range(phase.nSpecies()):
            self._spnames.append(phase.speciesName(k))
        _cantera.mix_addPhase(self.__mixid, phase.thermo_hndl(), moles)

    def nPhases(self):
        """Total number of phases defined for the mixture."""
        return len(self._phases)

    def phase(self, n):
        """Return the object representing the nth phase in the mixture."""
        return self._phases[n]

    def phaseName(self, n):
        """Name of phase *n*."""
        return self._phases[n].name()

    def phaseNames(self):
        """Names of all phases in the order added."""
        np = self.nPhases()
        nm = []
        for n in range(np):
            nm.append(self.phaseName(n))
        return nm

    def phaseIndex(self, phase):
        """Index of phase with name *phase*"""
        np = self.nPhases()
        if type(phase) <> types.StringType:
            return phase
        for n in range(np):
            if self.phaseName(n) == phase:
                return n
        return -1

    def nElements(self):
        """Total number of elements present in the mixture."""
        return _cantera.mix_nElements(self.__mixid)

    def elementIndex(self, element):
        """Index of element with name 'element'.

        >>> mix.elementIndex('H')
        2
        """
        if type(element) == types.StringType:
            return _cantera.mix_elementIndex(self.__mixid, element)
        else:
            return element

    def nSpecies(self):
        """Total number of species present in the mixture. This is the
        sum of the numbers of species in each phase."""
        return _cantera.mix_nSpecies(self.__mixid)

    def speciesName(self, k):
        """Name of the species with index *k*. Note that index numbers
        are assigned in order as phases are added."""
        return self._spnames[k]

    def speciesNames(self):
        n = self.nSpecies()
        s = []
        for k in range(n):
            s.append(self.speciesName(k))
        return s

    def speciesIndex(self, species):
        """Index of species with name *species*. If *species* is not a string,
        then it is simply returned."""
        if type(species) == types.StringType:
            return self._spnames.index(species)
        else:
            return species

    def nAtoms(self, k, m):
        """Number of atoms of element *m* in species *k*. Both the species and
        the element may be referenced either by name or by index number.

        >>> n = mix.nAtoms('CH4','H')
        4.0

        """
        kk = self.speciesIndex(k)
        mm = self.elementIndex(m)
        return _cantera.mix_nAtoms(self.__mixid, kk, mm)

    def setTemperature(self, t):
        """Set the temperature [K]. The temperatures of all phases are
        set to this value, holding the pressure fixed."""
        return _cantera.mix_setTemperature(self.__mixid, t)

    def temperature(self):
        """The temperature [K]."""
        return _cantera.mix_temperature(self.__mixid)

    def minTemp(self):
        """The minimum temperature for which all species in
        multi-species solutions have valid thermo data. Stoichiometric
        phases are not considered in determining minTemp.  """
        return _cantera.mix_minTemp(self.__mixid)

    def maxTemp(self):
        """The maximum temperature for which all species in
        multi-species solutions have valid thermo data. Stoichiometric
        phases are not considered in determining maxTemp.  """
        return _cantera.mix_maxTemp(self.__mixid)

    def charge(self):
        """The total charge in Coulombs, summed over all phases."""
        return _cantera.mix_charge(self.__mixid)

    def phaseCharge(self, p):
        """The charge of phase *p* (Coulombs)."""
        return _cantera.mix_phaseCharge(self.__mixid, p)

    def setPressure(self, p):
        """Set the pressure [Pa]. The pressures of all phases are set
        to the specified value, holding the temperature fixed."""
        return _cantera.mix_setPressure(self.__mixid, p)

    def pressure(self):
        """The pressure [Pa]."""
        return _cantera.mix_pressure(self.__mixid)

    def phaseMoles(self, n = -1):
        """Moles of phase *n*."""
        if n == -1:
            np = self.nPhases()
            moles = zeros(np,'d')
            for m in range(np):
                moles[m] = _cantera.mix_phaseMoles(self.__mixid, m)
            return moles
        else:
            return _cantera.mix_phaseMoles(self.__mixid, n)

    def setPhaseMoles(self, n, moles):
        """Set the number of moles of phase *n*."""
        _cantera.mix_setPhaseMoles(self.__mixid, n, moles)

    def setSpeciesMoles(self, moles):
        """Set the moles of the species [kmol]. The moles may be
        specified either as a string, or as an array. If an array is
        used, it must be dimensioned at least as large as the total
        number of species in the mixture. Note that the species may
        belong to any phase, and unspecified species are set to zero.

        >>> mix.setSpeciesMoles('C(s):1.0, CH4:2.0, O2:0.2')

        """
        if type(moles) == types.StringType:
            _cantera.mix_setMolesByName(self.__mixid, moles)
        else:
            _cantera.mix_setMoles(self.__mixid, asarray(moles))

    def speciesMoles(self, species = ""):
        """Moles of species k."""
        moles = zeros(self.nSpecies(),'d')
        for k in range(self.nSpecies()):
            moles[k] = _cantera.mix_speciesMoles(self.__mixid, k)
        return self.selectSpecies(moles, species)

    def elementMoles(self, m):
        """Total number of moles of element *m*, summed over all species.
        The element may be referenced either by index number or by name.
        """
        mm = self.elementIndex(m)
        return _cantera.mix_elementMoles(self.__mixid, mm)

    def chemPotentials(self, species=[]):
        """The chemical potentials of all species [J/kmol]."""
        mu = zeros(self.nSpecies(),'d')
        _cantera.mix_getChemPotentials(self.__mixid, mu)
        return self.selectSpecies(mu, species)

    def set(self, **p):
        for o in p.keys():
            v = p[o]
            if o == 'T' or o == 'Temperature':
                self.setTemperature(v)
            elif o == 'P' or o == 'Pressure':
                self.setPressure(v)
            elif o == 'Moles' or o == 'N':
                self.setSpeciesMoles(v)
            else:
                raise CanteraError("unknown property: "+o)

    def equilibrate(self, XY = "TP", err = 1.0e-9,
                    maxsteps = 1000, maxiter = 200, loglevel = 0):
        """Set the mixture to a state of chemical equilibrium.

        This method uses a version of the VCS algorithm to find the
        composition that minimizes the total Gibbs free energy of the
        mixture, subject to element conservation constraints. For a
        description of the theory, see Smith and Missen, "Chemical
        Reaction Equilibrium."  The VCS algorithm is implemented in
        Cantera kernel class ``MultiPhaseEquil``.

        The VCS algorithm solves for the equilibrium composition for
        specified temperature and pressure. If any other property pair
        other than ``TP`` is specified, then an outer iteration loop is
        used to adjust T and/or P so that the specified property
        values are obtained.

        :param XY:
            Two-letter string specifying the two properties to hold fixed.
            Currently, ``'TP'``, ``'HP'``, and ``'SP'`` are implemented.
            Default: ``'TP'``.
        :param err:
            Error tolerance. Iteration will continue until (Delta mu)/RT is
            less than this value for each reaction. Default: 1.0e-9. Note that
            this default is very conservative, and good equilibrium solutions
            may be obtained with larger error tolerances.
        :param maxsteps:
            Maximum number of steps to take while solving the equilibrium
            problem for specified *T* and *P*. Default: 1000.
        :param maxiter:
            Maximum number of temperature and/or pressure iterations.
            This is only relevant if a property pair other than (T,P) is
            specified. Default: 200.
        :param loglevel:
            Controls the amount of diagnostic output. If loglevel = 0, no
            diagnostic output is written. For values > 0, more detailed
            information is written to the log file as loglevel increases.
            The default is loglevel = 0.
            The logfile is written in HTML format, and may be viewed with
            any web browser. The default log file name is
            ``equilibrium_log.html``, but if this file exists, the log
            information will be written to "equilibrium_log{n}.html", where
            {n} is an integer chosen so that the log file does not already
            exist. Therefore, if 'equilibrate' is called multiple times,
            multiple log files will be written, with names
            ``equilibrate_log.html``, ``equilibrate_log1.html``,
            ``equilibrate_log2.html``, and so on. Existing log files will
            not be overwritten.

        >>> mix.equilibrate('TP')
        >>> mix.equilibrate('TP', err = 1.0e-6, maxiter = 500)

        """
        i = _cantera.mix_equilibrate(self.__mixid, XY, err, maxsteps,
                                        maxiter, loglevel)

    def vcs_equilibrate(self, XY = "TP", estimateEquil = 0, printLvl = 0,
                        solver = 2, rtol = 1.0e-9,
                        maxsteps = 1000, maxiter = 1000, loglevel = 0):
        """Set the mixture to a state of chemical equilibrium.

        This method uses a version of the VCS algorithm to find the
        composition that minimizes the total Gibbs free energy of the
        mixture, subject to element conservation constraints. For a
        description of the theory, see Smith and Missen, "Chemical
        Reaction Equilibrium."  The VCS algorithm is implemented in
        Cantera kernel class MultiPhaseEquil.

        The VCS algorithm solves for the equilibrium composition for
        specified temperature and pressure. If any other property pair
        other than ``'TP'`` is specified, then an outer iteration loop is
        used to adjust T and/or P so that the specified property
        values are obtained.

        :param XY:
            Two-letter string specifying the two properties to hold fixed.
            Currently, ``'TP'``, ``'HP'``, and ``'SP'`` are implemented.
            Default: ``'TP'``.
        :param printLvl:
            Controls the amount of diagnostic output written to cout. If
            printLvl = 0, no diagnostic output is written. For values > 0,
            more detailed information is written to cout.
            The default is printLvl = 0.
        :param solver:
            Determines which solver is used.
                - 1 MultiPhaseEquil solver
                - 2 VCSnonideal Solver (default)
        :param err:
            Error tolerance. Iteration will continue until (Delta mu)/RT is
            less than this value for each reaction. Default: 1.0e-9. Note that
            this default is very conservative, and good equilibrium solutions
            May be obtained with larger error tolerances.
        :param maxsteps:
            Maximum number of steps to take while solving the equilibrium
            problem for specified T and P. Default: 1000.
        :param maxiter:
            Maximum number of temperature and/or pressure iterations. This is
            only relevant if a property pair other than (T,P) is specified.
            Default: 200.
        :param loglevel:
            Controls the amount of diagnostic output written to html. If
            loglevel = 0, no diagnostic output is written. For values > 0,
            more detailed information is written to the log file as
            loglevel increases. The default is loglevel = 0.
            The logfile is written in HTML format, and may be viewed with
            any web browser. The default log file name is
            "equilibrium_log.html", but if this file exists, the log
            information will be written to "equilibrium_log{n}.html",
            where {n} is an integer chosen so that the log file does not
            already exist. Therefore, if 'equilibrate' is called multiple
            times, multiple log files will be written, with names
            "equilibrate_log.html", "equilibrate_log1.html",
            "equilibrate_log2.html", and so on. Existing log files will
            not be overwritten.
        """
        i = _cantera.mix_vcs_equilibrate(self.__mixid, XY, estimateEquil,
                                         printLvl, solver, rtol, maxsteps,
                                         maxiter, loglevel)

    def selectSpecies(self, f, species):
        """Given an array *f* of floating-point species properties,
        return an array of those values corresponding to species
        listed in *species*. This method is used internally to implement
        species selection in methods like :meth:`~.Phase.moleFractions`,
        :meth:`~.Phase.massFractions`, etc.

        >>> f = mix.chemPotentials()
        >>> muo2, muh2 = mix.selectSpecies(f, ['O2', 'H2'])
        """
        sp = []
        if species:
            if type(species) == types.StringType:
                sp = [species]
            else:
                sp = species
            fs = []
            k = 0
            for s in sp:
                k = self.speciesIndex(s)
                fs.append(f[k])
            return asarray(fs)
        else:
            return f
