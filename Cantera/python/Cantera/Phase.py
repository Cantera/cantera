"""Cantera.Phase

This module provides class Phase.

"""

import _cantera

import types
import Numeric 
from exceptions import CanteraError

__revision__ = "$Id$"

class Phase:

    """Class Phase manages basic state and constituent property
    information for a homogeneous phase of matter. It does not contain
    information on the equation of state, homogeneous kinetics, or
    transport properties -- these attributes of the phase are the
    responsibility of other classes (see the ThermoPhase, Kinetics, and
    Transport classes).

    In the C++ kernel, class Phase is only a base class for
    ThermoPhase. This structure is maintained in the Python wrapper
    class. Class Phase implements methods of the C++ Phase class, but
    the instance is actually a 'ThermoPhase' object. Also, in C++
    class Phase is itself derived from more elementary classes (State
    and Constituents), but these are not exposed in Python.

    """
    
    def __init__(self, index = -1):
        pass
        
    def phase_id(self):
        return self._phase_id
    
    def nElements(self):
        """Number of elements."""
        return _cantera.phase_nelements(self._phase_id)

    def atomicWeights(self, elements = []):
        """Array of element molar masses [kg/kmol]."""
        atw = _cantera.phase_getarray(self._phase_id,1)
        if elements:
            ae = []
            m = 0
            for e in elements:
                m = self.elementIndex(e)
                ae.append(atw[m])
            return Numeric.asarray(ae)
        else:
            return atw
        
    def nSpecies(self):
        """Number of species."""
        return _cantera.phase_nspecies(self._phase_id)

    def nAtoms(self, species = -1, element = -1):
        """Number of atoms of element 'element' in species 'species'.
        The  element and species may be specified by name or by number."""
        try:
            m = self.elementIndex(element)
            k = self.speciesIndex(species)
            na = _cantera.phase_natoms(self._phase_id, k, m)
            if na < 0: return 0
            return na
        except CanteraError:
            return 0
    
    def temperature(self):
        """Temperature [K]."""
        return _cantera.phase_temperature(self._phase_id)

    def density(self):
        """Mass density [kg/m^3]."""
        return _cantera.phase_density(self._phase_id)

    def volume_mass(self):
        """Specific volume [m^3/kg]."""
        return 1.0/_cantera.phase_density(self._phase_id)    

    def molarDensity(self):
        """Molar density [kmol/m^3]."""
        return _cantera.phase_molardensity(self._phase_id)

    def meanMolecularWeight(self):
        """Mean molar mass [kg/kmol].
        DEPRECATED: use meanMolarMass"""
        return _cantera.phase_meanmolwt(self._phase_id)

    def meanMolarMass(self):
        """Mean molar mass [kg/kmol]."""
        return _cantera.phase_meanmolwt(self._phase_id)

    def molarMasses(self, species = []):
        """Array of species molar masses [kg/kmol]."""
        mm  = _cantera.phase_getarray(self._phase_id,22)
        return self.selectSpecies(mm, species)

    def molecularWeights(self, species = []):
        """Array of species molar masses [kg/kmol].
        DEPRECATED: use molarMasses"""
        return self.molarMasses(species)

    def moleFractions(self, species = []):
        """Species mole fraction array."""
        x = _cantera.phase_getarray(self._phase_id,20)
        return self.selectSpecies(x, species)

    def moleFraction(self, species=-1):
        """Mole fraction of a species, referenced by name or
        index number.
        >>> ph.moleFraction(4)
        >>> ph.moleFraction('CH4')
        """
        k = self.speciesIndex(species)
        return _cantera.phase_molefraction(self._phase_id,k)    


    def massFractions(self, species = []):
        """Species mass fraction array."""
        y = _cantera.phase_getarray(self._phase_id,21)
        return self.selectSpecies(y, species)


    def massFraction(self, species=-1):
        """Mass fraction of one species, referenced by name or
        index number.
        >>> ph.massFraction(4)
        >>> ph.massFraction('CH4')
        """
        k = self.speciesIndex(species)
        return _cantera.phase_massfraction(self._phase_id,k)


    def elementName(self,m):
        """Name of element m."""
        return _cantera.phase_getstring(self._phase_id,1,m)

    def elementNames(self):
        """Return a tuple of all element names."""
        nel = self.nElements()
        return map(self.elementName,range(nel))

    def elementIndex(self, element=-1):
        """The index of element 'element', which may be specified as
        a string or an integer index. In the latter case, the index is
        checked for validity and returned. If no such element is
        present, an exception is thrown."""
        
        nel = self.nElements()
        if type(element) == types.IntType:
            m = element
        else:
            m = _cantera.phase_elementindex(self._phase_id, element)
        if m < 0 or m >= nel:
            raise CanteraError("""Element """+element+""" not in set """
                               +`self.elementNames()`)
        return m


    
    def speciesName(self,k):
        """Name of the species with index k."""
        return _cantera.phase_getstring(self._phase_id,2,k)

        
    def speciesNames(self):
        """Return a tuple of all species names."""        
        nsp = self.nSpecies()
        return map(self.speciesName,range(nsp))


    def speciesIndex(self, species=-1):
        """The index of species 'species', which may be specified as
        a string or an integer index. In the latter case, the index is
        checked for validity and returned. If no such species is
        present, an exception is thrown."""
        nsp = self.nSpecies()
        if type(species) == types.IntType or type(species) == types.FloatType:
            k = species
        else:
            k = _cantera.phase_speciesindex(self._phase_id,species)
        if k < 0 or k >= nsp:
            raise CanteraError("""Species """+`species`+""" not in set """
                               +`self.speciesNames()`)
        return k
        
    
    def setTemperature(self, t):
        """Set the temperature [K]."""
        _cantera.phase_setfp(self._phase_id,1,t)

    def setDensity(self, rho):
        """Set the density [kg/m3]."""
        _cantera.phase_setfp(self._phase_id,2,rho)
        
    def setMoleFractions(self, x, norm = 1):
        """Set the mole fractions. The values may be input either
        in a string or a sequence.
        >>> ph.setMoleFractions('CO:1, H2:7, H2O:7.8')
        >>> x = [1.0]*ph.nSpecies()
        >>> ph.setMoleFractions(x)
        By default, the input values will be scaled to sum to 1.0.
        If this is not desired, supply a third parameter 'norm' set to zero
        >>> ph.setMoleFractions(x, norm = 0)
        (Note that this only works if an array is input.)
        """
        if type(x) == types.StringType:
            _cantera.phase_setstring(self._phase_id,1,x)
        else:
            _cantera.phase_setarray(self._phase_id,1,norm,Numeric.asarray(x))

            
    def setMassFractions(self, x, norm = 1):
        """Set the mass fractions.
        See also: setMoleFractions
        """
        if type(x) == types.StringType:
            _cantera.phase_setstring(self._phase_id,2,x)
        else:        
            _cantera.phase_setarray(self._phase_id,2,norm,Numeric.asarray(x))        
        
    def setState_TRX(self, t, rho, x):
        """Set the temperature, density, and mole fractions."""
        self.setTemperature(t)
        self.setMoleFractions(x)
        self.setDensity(rho)

    def setState_TRY(self, t, rho, y):
        """Set the temperature, density, and mass fractions."""        
        self.setTemperature(t)
        self.setMassFractions(y)
        self.setDensity(rho)        

    def setState_TR(self, t, rho):
        """Set the temperature and density."""        
        self.setTemperature(t)
        self.setDensity(rho)
    
    def selectSpecies(self, f, sp):
        if sp:
            fs = []
            k = 0
            for s in sp:
                k = self.speciesIndex(s)
                fs.append(f[k])
            return Numeric.asarray(fs)
        else:
            return f
        
