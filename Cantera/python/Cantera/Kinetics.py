
from Cantera.exceptions import CanteraError, getCanteraError
from Cantera.ThermoPhase import ThermoPhase
from Cantera.XML import XML_Node
import Numeric
import _cantera


def buildKineticsPhases(root=None, id=None):
    """Return a list of ThermoPhase objects representing the phases
    involved in a reaction mechanism.

    root -- XML node contaning a 'kinetics' child
    id   -- id attribute of the desired 'kinetics' node
    """
    kin = root.child(id = id)
    phase_refs = kin.children("phaseRef")
    th = None
    phases = []
    for p in phase_refs:
        phase_id = p["id"]
        try:
            th = ThermoPhase(root=root, id=phase_id)
        except:
            if p["src"]:
                pnode = XML_Node(name="root",src=src)
                th = ThermoPhase(pnode, phase_id)
            else:
                raise CanteraError("phase "+phase_id+" not found.")
        phases.append(th)
    return phases


class Kinetics:
    """
    Kinetics managers. Instances of class Kinetics are responsible for
    evaluating reaction rates of progress, species production rates,
    and other quantities pertaining to a reaction mechanism.
    """
    
    def __init__(self, kintype=-1, thrm=0, xml_phase=None, id=None, phases=[]):
        """Build a kinetics manager from an XML specification.

        root     -- root of a CTML tree
        
        id       -- id of the 'kinetics' node within the tree that contains
                    the specification of the parameters.
        """
        np = len(phases)
        self._np = np
        self._ph = {}
        self._sp = []
        for p in phases:
            self._ph[p.thermophase()] = p
        self._phnum = {}
        
        self._end = [0]
        p0 = phases[0].thermophase()
        #self._ph[phases[0]] = phases
        #self._end.append(phases[0].nSpecies())
        p1 = -1
        p2 = -1
        p3 = -1
        p4 = -1
        if np >= 2:
            p1 = phases[1].thermophase()
        if np >= 3:
            p2 = phases[2].thermophase()
        if np >= 4:
            p3 = phases[3].thermophase()
        if np >= 5:
            p4 = phases[4].thermophase()
        if np >= 6:
            raise CanteraError("only 4 neighbor phases allowed")
        #try:
        self.ckin = _cantera.KineticsFromXML(xml_phase,
                                                 p0, p1, p2, p3, p4)
        
        for nn in range(self._np):
                p = self.phase(nn)
                self._phnum[p] = nn
                self._end.append(self._end[-1]+p.nSpecies())
                for k in range(p.nSpecies()):
                    self._sp.append(p.speciesName(k))
                
        #except:
        #    print getCanteraError()
        #    self.ckin = 0
        self.phases = phases
        

        
    def __del__(self):
        """Delete the kinetics manager."""
        if self.ckin > 0:
            _cantera.kin_delete(self.ckin)

    def kin_index(self):
        return self.ckin

    def kineticsType(self):
        """Kinetics manager type."""
        return _cantera.kin_type(self.ckin)

    def kineticsSpeciesIndex(self, name, phase):
        return _cantera.kin_speciesIndex(self.ckin, name, phase)

    def kineticsStart(self, n):
        return _cantera.kin_start(self.ckin, n)

    def phase(self, n):
        return self._ph[_cantera.kin_phase(self.ckin, n)]
    
    def nReactions(self):
        """Number of reactions."""
        return _cantera.kin_nreactions(self.ckin)

    def isReversible(self,i):
        """
        True (1) if reaction number 'i' is reversible,
        and false (0) otherwise.
        """
        return _cantera.kin_isreversible(self.ckin,i)

    def reactionType(self,i):
        """Type of reaction 'i'"""
        return _cantera.kin_rxntype(self.ckin,i)

    #def reactionString(self,i):
    #    return _cantera.kin_getstring(self.ckin,1,i)

    def reactionEqn(self,i):
        try:
            eqs = []
            for rxn in i:
                eqs.append(self.reactionString(rxn))
            return eqs
        except:
            return self.reactionString(i)

    def reactionString(self, i):
        s = ''
        nsp = _cantera.kin_nspecies(self.ckin)
        for k in range(nsp):
            nur = _cantera.kin_rstoichcoeff(self.ckin,k,i)
            if nur <> 0.0:
                if nur <> 1.0:
                    s += `int(nur)`+' '
                s += self._sp[k]+' + '
        s = s[:-2]
        if self.isReversible(i):
            s += ' <=> '
        else:
            s += ' => '
        for k in range(nsp):
            nup = _cantera.kin_pstoichcoeff(self.ckin,k,i)
            if nup <> 0.0:
                if nup <> 1.0:
                    s += `int(nup)`+' '
                s += self._sp[k]+' + '
        s = s[:-2]
        return s
    
    def reactantStoichCoeff(self,k,i):
        return _cantera.kin_rstoichcoeff(self.ckin,k,i)
    
    def reactantStoichCoeffs(self):
        nsp = _cantera.kin_nspecies(self.ckin)
        nr = _cantera.kin_nreactions(self.ckin)
        nu = Numeric.zeros((nsp,nr),'d')
        for i in range(nr):
            for k in range(nsp):
                nu[k,i] = _cantera.kin_rstoichcoeff(self.ckin,k,i)
        return nu

    def productStoichCoeff(self,k,i):
        return _cantera.kin_pstoichcoeff(self.ckin,k,i)    

    def productStoichCoeffs(self):
        nsp = _cantera.kin_nspecies(self.ckin)
        nr = _cantera.kin_nreactions(self.ckin)
        nu = Numeric.zeros((nsp,nr),'d')
        for i in range(nr):
            for k in range(nsp):
                nu[k,i] = _cantera.kin_pstoichcoeff(self.ckin,k,i)
        return nu
    
    def fwdRatesOfProgress(self):
        return _cantera.kin_getarray(self.ckin,10)

    def revRatesOfProgress(self):
        return _cantera.kin_getarray(self.ckin,20)

    def netRatesOfProgress(self):
        return _cantera.kin_getarray(self.ckin,30)

    def equilibriumConstants(self):
        return _cantera.kin_getarray(self.ckin,40)

    def creationRates(self, phase = None):
        c = _cantera.kin_getarray(self.ckin,50)
        if phase:
            if self._phnum.has_key(phase):
                n = self._phnum[phase]
                return c[self._end[n]:self._end[n+1]]
            else:
                raise CanteraError('unknown phase')
        else:
            return c                
                

    def destructionRates(self, phase = None):
        d = _cantera.kin_getarray(self.ckin,60)
        if phase:
            if self._phnum.has_key(phase):
                n = self._phnum[phase]
                return d[self._end[n]:self._end[n+1]]
            else:
                raise CanteraError('unknown phase')
        else:
            return d

    
    def netProductionRates(self, phase = None):
        w = _cantera.kin_getarray(self.ckin,70)        
        if phase:
            if self._phnum.has_key(phase):
                n = self._phnum[phase]
                return w[self._end[n]:self._end[n+1]]
            else:
                raise CanteraError('unknown phase')
        else:
            return w
        
    def sourceTerms(self):
        return _cantera.kin_getarray(self.ckin,80)        

    def multiplier(self,i):
        return _cantera.kin_multiplier(self.ckin,i)

    def setMultiplier(self,i,v):
        return _cantera.kin_setMultiplier(self.ckin,i,v)    
    
    def advanceCoverages(self,dt):
        return _cantera.kin_advanceCoverages(self.ckin,dt)    









    

    
