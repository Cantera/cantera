import _cantera
import types
from Numeric import zeros

class Mixture:
    """Class Mixture represents mixtures of one or more phases of matter."""
    
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
                    moles = 0
                self.addPhase(ph, moles)
                self._phases.append(ph)
                

    def __del__(self):
        _cantera.mix_del(self.__mixid)

    def __repr__(self):
        s = ''
        for p in range(len(self._phases)):
            s += '\n*******************    Phase '+`p`+'    ******************************\n'
            s += '\n Moles: '+`self.phaseMoles(p)`+'\n'
            s += self._phases[p].__repr__()+'\n\n'
        return s
    
    def addPhase(self, phase = None, moles = 0.0):
        for k in range(phase.nSpecies()):
            self._spnames.append(phase.speciesName(k))
        _cantera.mix_addPhase(self.__mixid, phase.thermo_hndl(), moles)
    
    def nElements(self):
        """Total number of elements present in the mixture."""
        return _cantera.mix_nElements(self.__mixid)
    def nSpecies(self):
        """Total number of species present in the mixture. This is the
        sum of the numbers of species in each phase."""
        return _cantera.mix_nSpecies(self.__mixid)
    def speciesName(self, k):
        return self._spnames[k]
    def speciesIndex(self, species):
        if type(species) == types.StringType:
            return self._spnames.index(species)
        else:
            return species
    def nAtoms(self, k, m):
        """Number of atoms of element m in species k."""
        return _cantera.mix_nAtoms(self.__mixid, k, m)
    def setTemperature(self, t):
        return _cantera.mix_setTemperature(self.__mixid, t)
    def temperature(self):
        return _cantera.mix_temperature(self.__mixid)
    def setPressure(self, p):
        return _cantera.mix_setPressure(self.__mixid, p)
    def pressure(self):
        return _cantera.mix_pressure(self.__mixid)    
    def phaseMoles(self, n):
        """Moles of phase n."""
        return _cantera.mix_phaseMoles(self.__mixid, n)
    def setPhaseMoles(self, n, moles):
        """Set the moles of phase n."""
        return _cantera.mix_setPhaseMoles(self.__mixid, n, moles)
    def speciesMoles(self, species):
        """Moles of species k."""
        k = self.speciesIndex(species)
        return _cantera.mix_speciesMoles(self.__mixid, k)
    def elementMoles(self, m):
        return _cantera.mix_elementMoles(self.__mixid, m)
    def chemPotentials(self):
        mu = zeros(self.nSpecies(),'d')
        _cantera.mix_getChemPotentials(self.__mixid, mu)
        return mu
    def equilibrate(self, XY = "TP", err = 1.0e-9, maxiter = 1000):
        return _cantera.mix_equilibrate(self.__mixid, XY, err, maxiter)
    
