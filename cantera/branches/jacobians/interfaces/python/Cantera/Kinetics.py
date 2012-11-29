"""
Kinetics managers.
"""

from Cantera.exceptions import CanteraError, getCanteraError
from Cantera.ThermoPhase import ThermoPhase
from Cantera.XML import XML_Node
from Cantera.num import zeros
import _cantera

class Kinetics:
    """
    Kinetics managers. Instances of class Kinetics are responsible for
    evaluating reaction rates of progress, species production rates,
    and other quantities pertaining to a reaction mechanism.
    """

    def __init__(self, kintype=-1, thrm=0, xml_phase=None, id=None, phases=[]):
        """
        Build a kinetics manager from an XML specification.

        :param kintype:
            Integer specifying the type of kinetics manager to create.
        :param root:
            Root of a CTML tree
        :param id:
            id of the 'kinetics' node within the tree that contains the
            specification of the parameters.
        """
        np = len(phases)
        self._sp = []
        self._phnum = {}

        # p0 through p4 are the integer indices of the phase objects
        # corresponding to the input sequence of phases
        self._end = [0]
        p0 = phases[0].thermophase()
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
            raise CanteraError("a maximum of 4 neighbor phases allowed")

        self.ckin = _cantera.KineticsFromXML(xml_phase,
                                                 p0, p1, p2, p3, p4)

        self._np = self.nPhases()
        for nn in range(self._np):
            p = self.phase(nn)
            self._phnum[p.thermophase()] = nn
            self._end.append(self._end[-1]+p.nSpecies())
            for k in range(p.nSpecies()):
                self._sp.append(p.speciesName(k))

    def __del__(self):
        self.clear()

    def clear(self):
        """Delete the kinetics manager."""
        if self.ckin > 0:
            _cantera.kin_delete(self.ckin)

    def kin_index(self):
        print "kin_index is deprecated. Use kinetics_hndl."
        return self.ckin

    def kinetics_hndl(self):
        return self.ckin

    def kineticsType(self):
        """Kinetics manager type."""
        return _cantera.kin_type(self.ckin)

    def kineticsSpeciesIndex(self, name, phase):
        """The index of a species.

        :param name:
            species name
        :param phase:
            phase name

        Kinetics managers for heterogeneous reaction mechanisms
        maintain a list of all species in all phases. The order of the
        species in this list determines the ordering of the arrays of
        production rates. This method returns the index for the
        specified species of the specified phase, and is used to
        locate the entry for a particular species in the production
        rate arrays.

        """
        return _cantera.kin_speciesIndex(self.ckin, name, phase)

    def kineticsStart(self, n):
        """The starting location of phase n in production rate arrays."""
        return _cantera.kin_start(self.ckin, n)

    def nPhases(self):
        """Number of phases."""
        return _cantera.kin_nPhases(self.ckin)

    def reactionPhaseIndex(self):
        """The phase in which the reactions take place."""
        return _cantera.kin_reactionPhaseIndex(self)

    def phase(self, n):
        """Return an object representing the nth phase."""
        return ThermoPhase(index = _cantera.kin_phase(self.ckin, n))

    def nReactions(self):
        """Number of reactions."""
        return _cantera.kin_nreactions(self.ckin)

    def isReversible(self,i):
        """
        True (1) if reaction number *i* is reversible,
        and false (0) otherwise.
        """
        return _cantera.kin_isreversible(self.ckin,i)

    def reactionType(self,i):
        """Type of reaction *i*"""
        return _cantera.kin_rxntype(self.ckin,i)

    def reactionEqn(self,i):
        """The equation for the specified reaction. If a list of equation numbers
        is given, then a list of equation strings is returned."""
        try:
            eqs = []
            for rxn in i:
                eqs.append(self.reactionString(rxn))
            return eqs
        except:
            return self.reactionString(i)

    def reactionString(self, i):
        """Reaction string for reaction number *i*"""
        s = ''
        nsp = _cantera.kin_nspecies(self.ckin)
        for k in range(nsp):
            nur = _cantera.kin_rstoichcoeff(self.ckin,k,i)
            if nur <> 0.0:
                if nur <> 1.0:
                    if nur <> round(nur):
                        s += str(nur)+' '
                    else:
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
                    if nup <> round(nup):
                        s += str(nup)+' '
                    else:
                        s += `int(nup)`+' '
                s += self._sp[k]+' + '
        s = s[:-2]
        return s

    def reactantStoichCoeff(self,k,i):
        """The stoichiometric coefficient of species *k* as a reactant in reaction *i*."""
        return _cantera.kin_rstoichcoeff(self.ckin,k,i)

    def reactantStoichCoeffs(self):
        """The array of reactant stoichiometric coefficients. Element
        [k,i] of this array is the reactant stoichiometric
        coefficient of species k in reaction i."""
        nsp = _cantera.kin_nspecies(self.ckin)
        nr = _cantera.kin_nreactions(self.ckin)
        nu = zeros((nsp,nr),'d')
        for i in range(nr):
            for k in range(nsp):
                nu[k,i] = _cantera.kin_rstoichcoeff(self.ckin,k,i)
        return nu

    def productStoichCoeff(self,k,i):
        """The stoichiometric coefficient of species *k* as a product in reaction *i*."""
        return _cantera.kin_pstoichcoeff(self.ckin,k,i)

    def productStoichCoeffs(self):
        """The array of product stoichiometric coefficients. Element
        [k,i] of this array is the product stoichiometric
        coefficient of species *k* in reaction *i*."""
        nsp = _cantera.kin_nspecies(self.ckin)
        nr = _cantera.kin_nreactions(self.ckin)
        nu = zeros((nsp,nr),'d')
        for i in range(nr):
            for k in range(nsp):
                nu[k,i] = _cantera.kin_pstoichcoeff(self.ckin,k,i)
        return nu

    def fwdRatesOfProgress(self):
        """Forward rates of progress of the reactions."""
        return _cantera.kin_getarray(self.ckin,10)

    def revRatesOfProgress(self):
        """Reverse rates of progress of the reactions."""
        return _cantera.kin_getarray(self.ckin,20)

    def netRatesOfProgress(self):
        """Net rates of progress of the reactions."""
        return _cantera.kin_getarray(self.ckin,30)

    def equilibriumConstants(self):
        """Equilibrium constants in concentration units for all reactions."""
        return _cantera.kin_getarray(self.ckin,40)

    def activationEnergies(self):
        """Activation energies in Kelvin for all reactions."""
        return _cantera.kin_getarray(self.ckin,32)

    def fwdRateConstants(self):
        """Forward rate constants for all reactions."""
        return _cantera.kin_getarray(self.ckin,34)

    def revRateConstants(self, doIrreversible = 0):
        """Reverse rate constants for all reactions."""
        if doIrreversible:
            return _cantera.kin_getarray(self.ckin,35)
        else:
            return _cantera.kin_getarray(self.ckin,36)

    def creationRates(self, phase = None):
        c = _cantera.kin_getarray(self.ckin,50)
        if phase:
            kp = phase.thermophase()
            if self._phnum.has_key(kp):
                n = self._phnum[kp]
                return c[self._end[n]:self._end[n+1]]
            else:
                raise CanteraError('unknown phase')
        else:
            return c


    def destructionRates(self, phase = None):
        d = _cantera.kin_getarray(self.ckin,60)
        if phase:
            kp = phase.thermophase()
            if self._phnum.has_key(kp):
                n = self._phnum[kp]
                return d[self._end[n]:self._end[n+1]]
            else:
                raise CanteraError('unknown phase')
        else:
            return d


    def netProductionRates(self, phase = None):
        w = _cantera.kin_getarray(self.ckin,70)
        if phase:
            kp = phase.thermophase()
            if self._phnum.has_key(kp):
                n = self._phnum[kp]
                return w[self._end[n]:self._end[n+1]]
            else:
                raise CanteraError('unknown phase')
        else:
            return w

    def sourceTerms(self):
        return _cantera.kin_getarray(self.ckin,80)

    def delta_H(self):
        return _cantera.kin_getarray(self.ckin,90)

    def delta_G(self):
        return _cantera.kin_getarray(self.ckin,91)

    def delta_S(self):
        return _cantera.kin_getarray(self.ckin,92)

    def delta_H0(self):
        return _cantera.kin_getarray(self.ckin,93)

    def delta_G0(self):
        return _cantera.kin_getarray(self.ckin,94)

    def delta_S0(self):
        return _cantera.kin_getarray(self.ckin,95)

    def multiplier(self,i):
        return _cantera.kin_multiplier(self.ckin,i)

    def setMultiplier(self, value = 0.0, reaction = -1):
        if reaction < 0:
            nr = self.nReactions()
            for i in range(nr):
                _cantera.kin_setMultiplier(self.ckin,i,value)
        else:
            _cantera.kin_setMultiplier(self.ckin,reaction,value)

    def advanceCoverages(self,dt):
        return _cantera.kin_advanceCoverages(self.ckin,dt)
