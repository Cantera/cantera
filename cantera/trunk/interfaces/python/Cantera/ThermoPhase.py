"""
 This module implements class ThermoPhase, a class representing
 thermodynamic phases.
"""
from Cantera.num import zeros
from Cantera.Phase import Phase

import _cantera
import types

class ThermoPhase(Phase):
    """
    A phase with an equation of state.

    Class ThermoPhase may be used to represent the intensive
    thermodynamic state of a phase of matter, which might be a gas,
    liquid, or solid.  Class ThermoPhase extends class Phase by
    providing methods that require knowledge of the equation of state.

    Class ThermoPhase is not usually instantiated directly. It is used
    as base class for classes :class:`~Cantera.Solution` and
    :class:`~Cantera.Interface.Interface`.
    """

    # used in the 'equilibrate' method
    _equilmap = {'TP':104,'TV':100,'HP':101,'SP':102,'SV':107,'UV':105,
                  'PT':104,'VT':100,'PH':101,'PS':102,'VS':107,'VU':105}


    def __init__(self, xml_phase=None, index=-1):
        """
        :param xml_phase:
            CTML node specifying the attributes of this phase
        :param index:
            optional. If positive, create only a Python wrapper for an existing
            kernel object, instead of creating a new kernel object. The value
            of *index* is the integer index number to reference the existing
            kernel object.
        """

        self._phase_id = 0
        self._owner = 0
        self.idtag = ""

        if index >= 0:
            # create a Python wrapper for an existing kernel
            # ThermoPhase instance
            self._phase_id = index

        elif xml_phase:
            # create a new kernel instance from an XML specification
            self._phase_id = _cantera.ThermoFromXML(xml_phase._xml_id)
            self.idtag = xml_phase["id"]
            self._owner = 1

        else:
            raise CanteraError('either xml_phase or index must be specified')


    def __del__(self):
        """Delete the object. If it is the owner of the kernel object,
        this is also deleted."""
        if self._owner:
            _cantera.thermo_delete(self._phase_id)

    def name(self):
        """The name assigned to the phase. The default value is the name
        attribute from the CTI file. But method setName can be used to
        set the name to anything desired, e.g. 'gas at inlet' or 'exhaust'
        """
        return self.idtag

    def setName(self, name):
        """ Set the name attribute. This can be any string"""
        self.idtag = name

    def refPressure(self):
        """Reference pressure [Pa].
        All standard-state thermodynamic properties are for this pressure.
        """
        return _cantera.thermo_refpressure(self._phase_id)

    def minTemp(self, sp=None):
        """ Minimum temperature for which thermodynamic property fits
        are valid.  If a species is specified (by name or number),
        then the minimum temperature is for only this
        species. Otherwise it is the lowest temperature for which the
        properties are valid for all species.  """
        if not sp:
            return _cantera.thermo_mintemp(self._phase_id, -1)
        else:
            return _cantera.thermo_mintemp(self._phase_id,
                                           self.speciesIndex(sp))

    def maxTemp(self, sp=None):
        """ Maximum temperature for which thermodynamic property fits
        are valid.  If a species is specified (by name or number),
        then the maximum temperature is for only this
        species. Otherwise it is the highest temperature for which the
        properties are valid for all species.  """
        if not sp:
            return _cantera.thermo_maxtemp(self._phase_id, -1)
        else:
            return _cantera.thermo_maxtemp(self._phase_id,
                                           self.speciesIndex(sp))

    def enthalpy_mole(self):
        """ The molar enthalpy [J/kmol]."""
        return _cantera.thermo_getfp(self._phase_id,1)

    def intEnergy_mole(self):
        """ The molar internal energy [J/kmol]."""
        return _cantera.thermo_getfp(self._phase_id,2)

    def entropy_mole(self):
        """ The molar entropy [J/kmol/K]."""
        return _cantera.thermo_getfp(self._phase_id,3)

    def gibbs_mole(self):
        """ The molar Gibbs function [J/kmol]."""
        return _cantera.thermo_getfp(self._phase_id,4)

    def cp_mole(self):
        """ The molar heat capacity at constant pressure [J/kmol/K]."""
        return _cantera.thermo_getfp(self._phase_id,5)

    def cv_mole(self):
        """ The molar heat capacity at constant volume [J/kmol/K]."""
        return _cantera.thermo_getfp(self._phase_id,6)

    def pressure(self):
        """ The pressure [Pa]."""
        return _cantera.thermo_getfp(self._phase_id,7)

    def electricPotential(self):
        """Electric potential [V]."""
        return _cantera.thermo_getfp(self._phase_id,25)

    def chemPotentials(self, species = []):
        """Species chemical potentials.

        This method returns an array containing the species
        chemical potentials [J/kmol]. The expressions used to
        compute these depend on the model implemented by the
        underlying kernel thermo manager."""
        mu = _cantera.thermo_getarray(self._phase_id,20)
        return self.selectSpecies(mu, species)

    def elementPotentials(self, elements = []):
        """Element potentials of the elements.

        This method returns an array containing the element potentials
        [J/kmol]. The element potentials are only defined for
        equilibrium states. This method first sets the composition to
        a state of equilibrium holding T and P constant, then computes
        the element potentials for this equilibrium state.  """

        lamb = _cantera.thermo_getarray(self._phase_id,21)
        return self.selectElements(lamb, elements)

    def enthalpies_RT(self, species = []):
        """Pure species non-dimensional reference state enthalpies.

        This method returns an array containing the pure-species
        standard-state enthalpies divided by RT. For gaseous species,
        these values are ideal gas enthalpies."""
        hrt = _cantera.thermo_getarray(self._phase_id,23)
        return self.selectSpecies(hrt, species)

    def entropies_R(self, species = []):
        """Pure species non-dimensional entropies.

        This method returns an array containing the pure-species
        standard-state entropies divided by R. For gaseous species,
        these values are ideal gas entropies."""
        sr = _cantera.thermo_getarray(self._phase_id,24)
        return self.selectSpecies(sr, species)

    def gibbs_RT(self, species = []):
        """Pure species non-dimensional Gibbs free energies.

        This method returns an array containing the pure-species
        standard-state Gibbs free energies divided by R.
        For gaseous species, these are ideal gas values."""
        grt = (_cantera.thermo_getarray(self._phase_id,23)
                - _cantera.thermo_getarray(self._phase_id,24))
        return self.selectSpecies(grt, species)

    def cp_R(self, species = []):
        """Pure species non-dimensional heat capacities
        at constant pressure.

        This method returns an array containing the pure-species
        standard-state heat capacities divided by R. For gaseous
        species, these values are ideal gas heat capacities."""
        cpr = _cantera.thermo_getarray(self._phase_id,25)
        return self.selectSpecies(cpr, species)


    def setPressure(self, p):
        """Set the pressure [Pa]."""
        _cantera.thermo_setfp(self._phase_id,1,p,0.0)

    def enthalpy_mass(self):
        """Specific enthalpy [J/kg]."""
        return _cantera.thermo_getfp(self._phase_id,8)

    def intEnergy_mass(self):
        """Specific internal energy [J/kg]."""
        return _cantera.thermo_getfp(self._phase_id,9)

    def entropy_mass(self):
        """Specific entropy [J/kg/K]."""
        return _cantera.thermo_getfp(self._phase_id,10)

    def gibbs_mass(self):
        """Specific Gibbs free energy [J/kg]."""
        return _cantera.thermo_getfp(self._phase_id,11)

    def cp_mass(self):
        """Specific heat at constant pressure [J/kg/K]."""
        return _cantera.thermo_getfp(self._phase_id,12)

    def cv_mass(self):
        """Specific heat at constant volume [J/kg/K]."""
        return _cantera.thermo_getfp(self._phase_id,13)

    def setState_TPX(self, t, p, x):
        """Set the temperature [K], pressure [Pa], and
        mole fractions."""
        self.setTemperature(t)
        self.setMoleFractions(x)
        self.setPressure(p)

    def setState_TPY(self, t, p, y):
        """Set the temperature [K], pressure [Pa], and
        mass fractions."""
        self.setTemperature(t)
        self.setMassFractions(y)
        self.setPressure(p)

    def setState_TP(self, t, p):
        """Set the temperature [K] and pressure [Pa]."""
        self.setTemperature(t)
        self.setPressure(p)

    def setState_PX(self, p, x):
        """Set the pressure [Pa], and mole fractions."""
        self.setMoleFractions(x)
        self.setPressure(p)

    def setState_PY(self, p, y):
        """Set the pressure [Pa], and mass fractions."""
        self.setMassFractions(y)
        self.setPressure(p)

    def setState_HP(self, h, p):
        """Set the state by specifying the specific enthalpy and
        the pressure."""
        _cantera.thermo_setfp(self._phase_id, 2, h, p)

    def setState_UV(self, u, v):
        """Set the state by specifying the specific internal
        energy and the specific volume."""
        _cantera.thermo_setfp(self._phase_id, 3, u, v)

    def setState_SV(self, s, v):
        """Set the state by specifying the specific entropy
        and the specific volume."""
        _cantera.thermo_setfp(self._phase_id, 4, s, v)

    def setState_SP(self, s, p):
        """Set the state by specifying the specific entropy
        energy and the pressure."""
        _cantera.thermo_setfp(self._phase_id, 5, s, p)

    def setElectricPotential(self, v):
        """Set the electric potential."""
        _cantera.thermo_setfp(self._phase_id, 6, v, 0);

    def equilibrate(self, XY, solver = -1, rtol = 1.0e-9,
                    maxsteps = 1000, maxiter = 100, loglevel = 0):
        """
        Set to a state of chemical equilibrium holding property pair
        *XY* constant.

        :param XY:
            A two-letter string, which must be one of the set::

                ['TP','TV','HP','SP','SV','UV','PT','VT','PH','PS','VS','VU']

            If H, U, S, or V is specified, the value must be the specific
            value (per unit mass)
        :param solver:
            Specifies the equilibrium solver to use. If solver = 0, a fast
            solver using the element potential method will be used. If
            solver > 0, a slower but more robust Gibbs minimization solver
            will be used. If solver < 0 or unspecified, the fast solver will
            be tried first, then if it fails the other will be tried.
        :param rtol:
            the relative error tolerance.
        :param maxsteps:
            maximum number of steps in composition to take to find a converged
            solution.
        :param maxiter:
            For the Gibbs minimization solver only, this specifies the number
            of 'outer' iterations on T or P when some property pair other than
            TP is specified.
        :param loglevel:
            Set to a value > 0 to write diagnostic output to a file in HTML
            format. Larger values generate more detailed information. The file
            will be named ``equilibrate_log.html.`` Subsequent files will be
            named ``equilibrate_log1.html``, etc., so that log files are
            not overwritten.
            """
        _cantera.thermo_equil(self._phase_id, XY, solver,
                              rtol, maxsteps, maxiter, loglevel)

    def saveState(self):
        """Return an array with state information that can later be
        used to restore the state."""
        state = zeros(self.nSpecies()+2,'d')
        state[0] = self.temperature()
        state[1] = self.density()
        state[2:] = self.massFractions()
        return state

    def restoreState(self, s):
        """Restore the state to that stored in array s."""
        self.setState_TRY(s[0], s[1], s[2:])

    def thermophase(self):
        """Return the integer index that is used to
        reference the kernel object. For internal use."""
        return self._phase_id

    def thermo_hndl(self):
        """Return the integer index that is used to
        reference the kernel object. For internal use."""
        return self._phase_id
