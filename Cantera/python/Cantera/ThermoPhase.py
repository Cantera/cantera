""" This module implements class ThermoPhase, a class representing
thermodynamic phases.  """

from Cantera.Phase import Phase

import _cantera
import types

def thermoIndex(id):
    return _cantera.thermo_thermoIndex(id)

class ThermoPhase(Phase):

    """ Class ThermoPhase may be used to represent the intensive state
    of a homogeneous phase of matter, which might be a gas, liquid, or solid.
    """

    #used in the 'equilibrate' method
    _equilmap = {'TP':104,'TV':100,'HP':101,'SP':102,'SV':107,'UV':105,
                  'PT':104,'VT':100,'PH':101,'PS':102,'VS':107,'VU':105}


    def __init__(self, xml_phase=None, index=-1):
        """Create a new object representing a phase of matter, or wrap
        an existing kernel instance."""
        
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

    def chemPotentials(self):
        """Species chemical potentials.
        
        This method returns an array containing the species
        chemical potentials [J/kmol]. The expressions used to
        compute these depend on the model implemented by the
        underlying kernel thermo manager."""
        return _cantera.thermo_getarray(self._phase_id,20)
    
    def enthalpies_RT(self):
        """Pure species non-dimensional enthalpies.
        
        This method returns an array containing the pure-species
        standard-state enthalpies divided by RT. For gaseous species,
        these values are ideal gas enthalpies."""
        return _cantera.thermo_getarray(self._phase_id,23)

    def entropies_R(self):
        """Pure species non-dimensional entropies.
        
        This method returns an array containing the pure-species
        standard-state entropies divided by R. For gaseous species,
        these values are ideal gas entropies."""
        return _cantera.thermo_getarray(self._phase_id,24)    

    def gibbs_RT(self):
        """Pure species non-dimensional Gibbs free energies.
        
        This method returns an array containing the pure-species
        standard-state Gibbs free energies divided by R.
        For gaseous species, these are ideal gas values."""        
        return (_cantera.thermo_getarray(self._phase_id,23)
                - _cantera.thermo_getarray(self._phase_id,24))
    
    def cp_R(self):
        """Pure species non-dimensional heat capacities
        at constant pressure.
        
        This method returns an array containing the pure-species
        standard-state heat capacities divided by R. For gaseous
        species, these values are ideal gas heat capacities."""
        return _cantera.thermo_getarray(self._phase_id,25)
    

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
        _cantera.thermo_setfp(self._phase_id, 6, v, 0);
        
    def equilibrate(self, XY):
        """Set to a state of chemical equilibrium holding property pair
        'XY' constant. The pair is specified by a two-letter string,
        which must be one of the set
        ['TP','TV','HP','SP','SV','UV','PT','VT','PH','PS','VS','VU'].
        If H, U, S, or V is specified, the value must be the specific
        value (per unit mass).
        """
        ixy = ThermoPhase._equilmap[XY]
        if ixy > 0:
            _cantera.thermo_equil(self._phase_id, ixy)
        else:
            raise 'invalid equilibrium option: '+XY
        

    def critTemperature(self):
        """Critical temperature [K]."""
        return _cantera.thermo_getfp(self._phase_id,50)

    def critPressure(self):
        """Critical pressure [Pa]."""
        return _cantera.thermo_getfp(self._phase_id,51)

    def critDensity(self):
        """Critical density [kg/m3]."""
        return _cantera.thermo_getfp(self._phase_id,52)

    def vaporFraction(self):
        """Vapor fraction."""
        return _cantera.thermo_getfp(self._phase_id,53)

    def setState_Psat(self, p, vaporFraction):
        _cantera.thermo_setfp(self._phase_id,8, p, vaporFraction)

    def setState_Tsat(self, t, vaporFraction):
        _cantera.thermo_setfp(self._phase_id,7, t, vaporFraction)        
        

    def thermophase(self):
        """Return the integer index that is used to
        reference the kernel object. For internal use."""
        return self._phase_id

    def thermo_hndl(self):
        """Return the integer index that is used to
        reference the kernel object. For internal use."""        
        return self._phase_id
    



    

    
