# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import warnings

# Need a pure-python class to store weakrefs to
class _WeakrefProxy:
    pass

cdef class Mixture:
    """

    Class Mixture represents mixtures of one or more phases of matter.  To
    construct a mixture, supply a list of phases to the constructor, each
    paired with the number of moles for that phase::

        >>> gas = cantera.Solution('gas.cti')
        >>> gas.species_names
        ['H2', 'H', 'O2', 'O', 'OH']
        >>> graphite = cantera.Solution('graphite.cti')
        >>> graphite.species_names
        ['C(g)']
        >>> mix = cantera.Mixture([(gas, 1.0), (graphite, 0.1)])
        >>> mix.species_names
        ['H2', 'H', 'O2', 'O', 'OH', 'C(g)']

    Note that the objects representing each phase compute only the intensive
    state of the phase -- they do not store any information on the amount of
    this phase. Mixture objects, on the other hand, represent the full
    extensive state.

    Mixture objects are 'lightweight' in the sense that they do not store
    parameters needed to compute thermodynamic or kinetic properties of the
    phases. These are contained in the ('heavyweight') phase objects. Multiple
    mixture objects may be constructed using the same set of phase objects.
    Each one stores its own state information locally, and synchronizes the
    phases objects whenever it requires phase properties.
    """
    def __cinit__(self, phases):
        self.mix = new CxxMultiPhase()
        self._phases = []
        self._weakref_proxy = _WeakrefProxy()

        cdef _SolutionBase phase
        if isinstance(phases[0], _SolutionBase):
            # Assign default composition to the list of phases
            phases = [(p, 1 if i == 0 else 0) for i,p in enumerate(phases)]

        for phase,moles in phases:
            phase._references[self._weakref_proxy] = True
            self.mix.addPhase(phase.thermo, moles)
            self._phases.append(phase)

        self.mix.init()
        if self._phases:
            self.P = self._phases[0].P
            self.T = self._phases[0].T

    def __dealloc__(self):
        del self.mix

    def report(self, threshold=1e-14):
        """
        Generate a report describing the thermodynamic state of this mixture. To
        print the report to the screen, simply call the mixture object. The
        following two statements are equivalent::

        >>> mix()
        >>> print(mix.report())
        """
        self.mix.updatePhases()
        s = []
        for i,phase in enumerate(self._phases):
            s.append('************ Phase {0} ************'.format(phase.name))
            s.append('Moles:  {0}'.format(self.phase_moles(i)))
            s.append(phase.report(threshold=threshold))

        return '\n'.join(s)

    def __call__(self):
        print(self.report())

    property n_elements:
        """Total number of elements present in the mixture."""
        def __get__(self):
            return self.mix.nElements()

    cpdef int element_index(self, element) except *:
        """Index of element with name 'element'.

            >>> mix.element_index('H')
            2
        """
        if isinstance(element, (str, bytes)):
            index = self.mix.elementIndex(stringify(element))
        elif isinstance(element, (int, float)):
            index = <int>element
        else:
            raise TypeError("'element' must be a string or a number")

        if not 0 <= index < self.n_elements:
            raise ValueError('No such element.')

        return index

    property n_species:
        """Number of species."""
        def __get__(self):
            return self.mix.nSpecies()

    def species_name(self, k):
        """Name of the species with index *k*. Note that index numbers
        are assigned in order as phases are added."""
        return pystr(self.mix.speciesName(k))

    property species_names:
        def __get__(self):
            return [self.species_name(k) for k in range(self.n_species)]

    def species_index(self, phase, species):
        """
        :param phase:
            Phase object, index or name
        :param species:
            Species name or index

        Returns the global index of species *species* in phase *phase*.
        """
        p = self.phase_index(phase)

        if isinstance(species, (str, bytes)):
            k = self.phase(p).species_index(species)
        elif isinstance(species, (int, float)):
            k = <int?>species
            if not 0 <= k < self.n_species:
                raise ValueError('Species index out of range')
        else:
            raise TypeError("'species' must be a string or number")

        return self.mix.speciesIndex(k, p)

    def n_atoms(self, k, m):
        """
        Number of atoms of element *m* in the species with global index *k*.
        The element may be referenced either by name or by index.

        >>> n = mix.n_atoms(3, 'H')
        4.0
        """
        if not 0 <= k < self.n_species:
            raise IndexError('Species index ({0}) out of range (0 < {1})'.format(k, self.n_species))
        return self.mix.nAtoms(k, self.element_index(m))

    property n_phases:
        """Number of phases"""
        def __get__(self):
            return len(self._phases)

    def phase(self, n):
        return self._phases[n]

    def phase_index(self, p):
        """Index of the phase named *p*."""
        if isinstance(p, ThermoPhase):
            p = p.name

        if isinstance(p, (int, float)):
            if p == int(p) and 0 <= p < self.n_phases:
                return int(p)
            else:
                raise IndexError("Phase index '{0}' out of range.".format(p))
        elif isinstance(p, (str, bytes)):
            for i, phase in enumerate(self._phases):
                if phase.name == p:
                    return i
        raise KeyError("No such phase: '{0}'".format(p))

    property phase_names:
        """Names of all phases in the order added."""
        def __get__(self):
            return [phase.name for phase in self._phases]

    property T:
        """
        Get or set the Temperature [K] of all phases in the mixture. When set,
        the pressure of the mixture is held fixed.
        """
        def __get__(self):
            return self.mix.temperature()
        def __set__(self, T):
            self.mix.setTemperature(T)

    property min_temp:
        """
        The minimum temperature for which all species in multi-species
        solutions have valid thermo data. Stoichiometric phases are not
        considered in determining min_temp.
        """
        def __get__(self):
            return self.mix.minTemp()

    property max_temp:
        """
        The maximum temperature for which all species in multi-species
        solutions have valid thermo data. Stoichiometric phases are not
        considered in determining max_temp.
        """
        def __get__(self):
            return self.mix.maxTemp()

    property P:
        """Get or set the Pressure [Pa] of all phases in the mixture. When set,
         the temperature of the mixture is held fixed."""
        def __get__(self):
            return self.mix.pressure()
        def __set__(self, P):
            self.mix.setPressure(P)

    property charge:
        """The total charge in Coulombs, summed over all phases."""
        def __get__(self):
            return self.mix.charge()

    def phase_charge(self, p):
        """The charge of phase *p* in Coulombs."""
        return self.mix.phaseCharge(self.phase_index(p))

    def phase_moles(self, p=None):
        """
        Moles in phase *p*, if *p* is specified, otherwise the number of
        moles in all phases.
        """
        if p is None:
            return [self.mix.phaseMoles(n) for n in range(self.n_phases)]
        else:
            return self.mix.phaseMoles(self.phase_index(p))

    def set_phase_moles(self, p, moles):
        """
        Set the number of moles of phase *p* to *moles*
        """
        self.mix.setPhaseMoles(self.phase_index(p), moles)

    property species_moles:
        """
        Get or set the number of moles of each species. May be set either as a
        string or as an array. If an array is used, it must be dimensioned at
        least as large as the total number of species in the mixture. Note that
        the species may belong to any phase, and unspecified species are set to
        zero.

        >>> mix.species_moles = 'C(s):1.0, CH4:2.0, O2:0.2'
        """
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(self.n_species)
            for k in range(self.n_species):
                data[k] = self.mix.speciesMoles(k)
            return data

        def __set__(self, moles):
            if isinstance(moles, (str, bytes)):
                self.mix.setMolesByName(stringify(moles))
                return

            if len(moles) != self.n_species:
                raise ValueError('mole array must be of length n_species')

            cdef np.ndarray[np.double_t, ndim=1] data = \
                np.ascontiguousarray(moles, dtype=np.double)
            self.mix.setMoles(&data[0])

    def element_moles(self, e):
        """
        Total number of moles of element *e*, summed over all species.
        The element may be referenced either by index number or by name.
        """
        return self.mix.elementMoles(self.element_index(e))

    property chemical_potentials:
        """The chemical potentials of all species [J/kmol]."""
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(self.n_species)
            self.mix.getChemPotentials(&data[0])
            return data

    def equilibrate(self, XY, solver='auto', rtol=1e-9, max_steps=1000,
                    max_iter=100, estimate_equil=0, log_level=0):
        """
        Set to a state of chemical equilibrium holding property pair *XY*
        constant. This method uses a version of the VCS algorithm to find the
        composition that minimizes the total Gibbs free energy of the mixture,
        subject to element conservation constraints. For a description of the
        theory, see Smith and Missen, "Chemical Reaction Equilibrium."

        :param XY:
            A two-letter string, which must be one of the set::

                ['TP', 'HP', 'SP']
        :param solver: Set to either 'auto', 'vcs', or 'gibbs' to choose
            implementation of the solver to use. 'vcs' uses the solver
            implemented in the C++ class 'VCSnonideal', 'gibbs' uses the one
            implemented in class 'MultiPhaseEquil'. 'auto' will try the 'vcs'
            solver first and then the 'gibbs' solver if that fails.
        :param rtol:
            Error tolerance. Iteration will continue until (Delta mu)/RT is
            less than this value for each reaction. Note that this default is
            very conservative, and good equilibrium solutions may be obtained
            with larger error tolerances.
        :param max_steps:
            Maximum number of steps to take while solving the equilibrium
            problem for specified *T* and *P*.
        :param max_iter:
            Maximum number of temperature and/or pressure iterations.
            This is only relevant if a property pair other than (T,P) is
            specified.
        :param estimate_equil:
            Flag indicating whether the solver should estimate its own initial
            condition. If 0, the initial mole fraction vector in the phase
            objects are used as the initial condition. If 1, the initial mole
            fraction vector is used if the element abundances are satisfied.
            if -1, the initial mole fraction vector is thrown out, and an
            estimate is formulated.
        :param log_level:
            Determines the amount of output displayed during the solution
            process. 0 indicates no output, while larger numbers produce
            successively more verbose information.
        """
        self.mix.equilibrate(stringify(XY.upper()), stringify(solver), rtol,
                             max_steps, max_iter, estimate_equil, log_level)
