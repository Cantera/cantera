"""
Module Thermo
"""

DEPRECATED 


from Cantera.Phase import Phase

import ctthermo
import ctphase
import types

def thermoIndex(id):
    return ctthermo.thermoIndex(id)

class Thermo:

    _equilmap = {'TP':104,'TV':100,'HP':101,'SP':102,'SV':107,'UV':105,
                  'PT':104,'VT':100,'PH':101,'PS':102,'VS':107,'VU':105}

    def __init__(self, eostype=1, phase=None, sptherm=0,
                 root=None, id=None, index=-1):
        self.__phase = None
        self.cthermo = None
        self._owner = 1
        self.idtag = ""

        #if thermoIndex(id) > 0:
        #    index = thermoIndex(id)
            
        if index >= 0:
            # create a Python wrapper for an existing kernel
            # Thermo instance 
            self.cthermo = index
            self.__phase = Phase(index = ctthermo.phase(index))
            self._owner = 0
                                 
        elif root:
            # create a new kernel instance from an XML specification
            self.cthermo, ph = ctthermo.ThermoFromXML(root._xml_id, id)
            self.__phase = Phase(ph)
            self.idtag = id
            
        else:
            # create a new kernel instance with specified parameters
            self.cthermo = ctthermo.Thermo(eostype, phase.phase_id(), sptherm)
            self.__phase = phase

    def __del__(self):
        if self._owner:
            ctthermo.delete(self.cthermo)

    def importFromXML(self, xml_root, id):
        ctthermo.import_xml(self.cthermo, xml_root._xml_id, id)
        
    def thermophase(self):
        return self.__phase

    def refPressure(self):
        """Reference pressure [Pa].
        All standard-state thermodynamic properties are for this pressure.
        """
        return ctthermo.refpressure(self.cthermo)

    def minTemp(self, sp=-1):
        """ Minimum temperature for which the parameterization of
        standard-state thermodynamic properties vs. T for species 'sp'
        is valid. If no species is specified, the value returned is
        the maximum value of minTemp for any one species, and
        therefore is the minimum temperature at which mixture
        thermodynamic properties are valid."""
        return ctthermo.mintemp(self.cthermo, self.speciesIndex(sp))

    def maxTemp(self, sp=-1):
        """ Maximum temperature for which the parameterization of
        standard-state thermodynamic properties vs. T for species 'sp'
        is valid. If no species is specified, the value returned is
        the minimum value of maxTemp for any one species, and
        therefore is the maximum temperature at which mixture
        thermodynamic properties are valid."""        
        return ctthermo.maxtemp(self.cthermo, self.speciesIndex(sp))    

    def enthalpy_mole(self):
        """ The molar enthalpy [J/kmol]."""
        return ctthermo.getfp(self.cthermo,1)

    def intEnergy_mole(self):
        """ The molar internal energy [J/kmol]."""
        return ctthermo.getfp(self.cthermo,2)          

    def entropy_mole(self):
        """ The molar entropy [J/kmol/K]."""
        return ctthermo.getfp(self.cthermo,3)

    def gibbs_mole(self):
        """ The molar Gibbs function [J/kmol]."""
        return ctthermo.getfp(self.cthermo,4)
    
    def cp_mole(self):
        """ The molar heat capacity at constant pressure [J/kmol/K]."""
        return ctthermo.getfp(self.cthermo,5)

    def cv_mole(self):
        """ The molar heat capacity at constant volume [J/kmol/K]."""
        return ctthermo.getfp(self.cthermo,6)

    def pressure(self):
        """ The pressure [Pa]."""
        return ctthermo.getfp(self.cthermo,7)        

    def chemPotentials(self):
        """Species chemical potentials.
        
        This method returns an array containing the species
        chemical potentials [J/kmol]. The expressions used to
        compute these depend on the model implemented by the
        underlying kernel thermo manager."""
        return ctthermo.getarray(self.cthermo,20)
    
    def enthalpies_RT(self):
        """Pure species non-dimensional enthalpies.
        
        This method returns an array containing the pure-species
        standard-state enthalpies divided by RT. For gaseous species,
        these values are ideal gas enthalpies."""
        return ctthermo.getarray(self.cthermo,23)

    def entropies_R(self):
        """Pure species non-dimensional entropies.
        
        This method returns an array containing the pure-species
        standard-state entropies divided by R. For gaseous species,
        these values are ideal gas entropies."""
        return ctthermo.getarray(self.cthermo,24)    

    def gibbs_RT(self):
        """Pure species non-dimensional Gibbs free energies.
        
        This method returns an array containing the pure-species
        standard-state Gibbs free energies divided by R.
        For gaseous species, these are ideal gas values."""        
        return (ctthermo.getarray(self.cthermo,23)
                - ctthermo.getarray(self.cthermo,24))
    
    def cp_R(self):
        """Pure species non-dimensional heat capacities
        at constant pressure.
        
        This method returns an array containing the pure-species
        standard-state heat capacities divided by R. For gaseous
        species, these values are ideal gas heat capacities."""
        return ctthermo.getarray(self.cthermo,25)
    

    def setPressure(self, p):
        """Set the pressure [Pa]."""
        ctthermo.setfp(self.cthermo,1,p,0.0)

    def enthalpy_mass(self):
        """Specific enthalpy [J/kg]."""
        return ctthermo.getfp(self.cthermo,8)        

    def intEnergy_mass(self):
        """Specific internal energy [J/kg]."""
        return ctthermo.getfp(self.cthermo,9)        

    def entropy_mass(self):
        """Specific entropy [J/kg/K]."""
        return ctthermo.getfp(self.cthermo,10)        

    def gibbs_mass(self):
        """Specific Gibbs free energy [J/kg]."""
        return ctthermo.getfp(self.cthermo,11)        

    def cp_mass(self):
        """Specific heat at constant pressure [J/kg/K]."""
        return ctthermo.getfp(self.cthermo,12)        

    def cv_mass(self):
        """Specific heat at constant volume [J/kg/K]."""        
        return ctthermo.getfp(self.cthermo,13)        

    def setState_HP(self, h, p):
        """Set the state by specifying the specific enthalpy and
        the pressure."""
        ctthermo.setfp(self.cthermo, 2, h, p)

    def setState_UV(self, u, v):
        """Set the state by specifying the specific internal
        energy and the specific volume."""
        ctthermo.setfp(self.cthermo, 3, u, v)        

    def setState_SV(self, s, v):
        """Set the state by specifying the specific entropy
        and the specific volume."""        
        ctthermo.setfp(self.cthermo, 4, s, v)
        
    def setState_SP(self, s, p):
        """Set the state by specifying the specific entropy
        energy and the pressure."""        
        ctthermo.setfp(self.cthermo, 5, s, p)        

    def equilibrate(self, XY):
        """Set to a state of chemical equilibrium holding property pair
        'XY' constant. The pair is specified by a two-letter string,
        which must be one of the set
        ['TP','TV','HP','SP','SV','UV','PT','VT','PH','PS','VS','VU'].
        If H, U, S, or V is specified, the value must be the specific
        value (per unit mass).
        """
        ixy = Thermo._equilmap[XY]
        if ixy > 0:
            ctthermo.equil(self.cthermo, ixy)
        else:
            raise 'invalid equilibrium option: '+XY
        







    

    
