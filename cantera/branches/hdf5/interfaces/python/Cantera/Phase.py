"""Cantera.Phase

This module provides class Phase.

"""

import _cantera

import types
from Cantera.num import asarray
from exceptions import CanteraError

# return true is x is a sequence
def _isseq(n, x):
    try:
        y = x[n-1]
        return 1
    except:
        return 0

class Phase:

    """Phases of matter.

    Class Phase manages basic state and constituent property
    information for a homogeneous phase of matter. It handles only
    those properties that do not require the equation of state, namely
    the temperature, density, chemical composition, and attributes of
    the elements and species.

    It does not know about the pressure, or any other thermodynamic property
    requiring the equation of state -- class ThermoPhase derives from Phase
    and adds those properties.

    Class Phase is not usually instantiated directly. It is used as a
    base class for class ThermoPhase.

    """

    #def __init__(self, index = -1):
    #    pass

    def phase_id(self):
        """The integer index used to access the kernel-level object.
        Internal."""
        return self._phase_id

    def nElements(self):
        """Number of elements."""
        return _cantera.phase_nelements(self._phase_id)

    def atomicWeights(self, elements = []):
        """Array of element molar masses [kg/kmol].

        If a sequence of element symbols is supplied, only the values
        for those elements are returned, ordered as in the
        list. Otherwise, the values are for all elements in the phase,
        ordered as in the input file.  """
        atw = _cantera.phase_getarray(self._phase_id,1)
        if elements:
            ae = []
            m = 0
            for e in elements:
                m = self.elementIndex(e)
                ae.append(atw[m])
            return asarray(ae)
        else:
            return atw

    def nSpecies(self):
        """Number of species."""
        return _cantera.phase_nspecies(self._phase_id)

    def nAtoms(self, species = None, element = None):
        """Number of atoms of element *element* in species *species*.
        The  element and species may be specified by name or by number.

        >>> ph.nAtoms('CH4','H')
        4

        """
        try:
            m = self.elementIndex(element)
            k = self.speciesIndex(species)
            na = _cantera.phase_natoms(self._phase_id, k, m)
            #if na < 0: return 0
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
        """Mean molar mass [kg/kmol]."""
        return _cantera.phase_meanmolwt(self._phase_id)

    def meanMolarMass(self):
        """Mean molar mass [kg/kmol]."""
        return _cantera.phase_meanmolwt(self._phase_id)

    def molarMasses(self, species = None):
        """Array of species molar masses [kg/kmol]."""
        mm  = _cantera.phase_getarray(self._phase_id,22)
        return self.selectSpecies(mm, species)

    def molecularWeights(self, species = None):
        """Array of species molar masses [kg/kmol]."""
        return self.molarMasses(species)

    def moleFractions(self, species = None):
        """Species mole fraction array.
        If optional argument *species* is supplied, then only the values
        for the selected species are returned.

        >>> x1 = ph.moleFractions()   # all species
        >>> x2 = ph.moleFractions(['OH', 'CH3'. 'O2'])
        """
        x = _cantera.phase_getarray(self._phase_id,20)
        return self.selectSpecies(x, species)

    def moleFraction(self, species):
        """Mole fraction of a species, referenced by name or index number.

        >>> ph.moleFraction(4)
        >>> ph.moleFraction('CH4')
        """
        k = self.speciesIndex(species)
        return _cantera.phase_molefraction(self._phase_id,k)


    def massFractions(self, species = None):
        """Species mass fraction array.
        If optional argument *species* is supplied, then only the values for
        the selected species are returned.

        >>> y1 = ph.massFractions()   # all species
        >>> y2 = ph.massFractions(['OH', 'CH3'. 'O2'])
        """
        y = _cantera.phase_getarray(self._phase_id,21)
        return self.selectSpecies(y, species)


    def massFraction(self, species):
        """Mass fraction of one species, referenced by name or
        index number.

        >>> ph.massFraction(4)
        >>> ph.massFraction('CH4')
        """
        k = self.speciesIndex(species)
        return _cantera.phase_massfraction(self._phase_id,k)


    def elementName(self,m):
        """Name of the element with index number *m*."""
        return _cantera.phase_getstring(self._phase_id,1,m)

    def elementNames(self):
        """Return a tuple of all element names."""
        nel = self.nElements()
        return map(self.elementName,range(nel))

    def elementIndex(self, element):
        """The index of element *element*, which may be specified as
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
        """Name of the species with index *k*."""
        return _cantera.phase_getstring(self._phase_id,2,k)


    def speciesNames(self):
        """Return a tuple of all species names."""
        nsp = self.nSpecies()
        return map(self.speciesName,range(nsp))


    def speciesIndex(self, species):
        """The index of species *species*, which may be specified as
        a string or an integer index. In the latter case, the index is
        checked for validity and returned. If no such species is
        present, an exception is thrown."""
        nsp = self.nSpecies()
        if type(species) == types.ListType:
            s = []
            for sp in species:
                s.append(self.speciesIndex(sp))
            return s

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

    def setMolarDensity(self, n):
        """Set the density [kmol/m3]."""
        _cantera.phase_setfp(self._phase_id,3,n)

    def setMoleFractions(self, x, norm = 1):
        """Set the mole fractions.

        :param x:
            string or array of mole fraction values
        :param norm:
            If non-zero (default), array values will be scaled to sum to 1.0.

        >>> ph.setMoleFractions('CO:1, H2:7, H2O:7.8')
        >>> x = [1.0]*ph.nSpecies()
        >>> ph.setMoleFractions(x)
        >>> ph.setMoleFractions(x, norm = 0)  # don't normalize values
        """
        if type(x) == types.StringType:
            _cantera.phase_setstring(self._phase_id,1,x)
        elif _isseq(self.nSpecies(), x):
            _cantera.phase_setarray(self._phase_id,1,norm,asarray(x))
        else:
            raise CanteraError('mole fractions must be a string or array')


    def setMassFractions(self, x, norm = 1):
        """Set the mass fractions.
        See :meth:`~.Phase.setMoleFractions`
        """
        if type(x) == types.StringType:
            _cantera.phase_setstring(self._phase_id,2,x)
        elif _isseq(self.nSpecies(), x):
            _cantera.phase_setarray(self._phase_id,2,norm,asarray(x))
        else:
            raise CanteraError('mass fractions must be a string or array')


    def setState_TRX(self, t, rho, x):
        """Set the temperature, density, and mole fractions. The mole
        fractions may be entered as a string or array,
        >>> ph.setState_TRX(600.0, 2.0e-3, 'CH4:0.4, O2:0.6')
        """

        self.setTemperature(t)
        self.setMoleFractions(x)
        self.setDensity(rho)

    def setState_TNX(self, t, n, x):
        """Set the temperature, molardensity, and mole fractions. The mole
        fractions may be entered as a string or array,

        >>> ph.setState_TNX(600.0, 2.0e-3, 'CH4:0.4, O2:0.6')
        """

        self.setTemperature(t)
        self.setMoleFractions(x)
        self.setMolarDensity(n)

    def setState_TRY(self, t, rho, y):
        """Set the temperature, density, and mass fractions."""
        self.setTemperature(t)
        self.setMassFractions(y)
        self.setDensity(rho)

    def setState_TR(self, t, rho):
        """Set the temperature and density, leaving the composition
        unchanged."""
        self.setTemperature(t)
        self.setDensity(rho)

    def selectSpecies(self, f, species):
        """Given an array *f* of floating-point species properties, return
        those values corresponding to species listed in *species*. Returns an
        array if *species* is a sequence, or a scalar if *species* is a
        scalar. This method is used internally to implement species selection
        in methods like moleFractions, massFractions, etc.

        >>> f = ph.chemPotentials()
        >>> muo2, muh2 = ph.selectSpecies(f, ['O2', 'H2'])
        >>> muh2 = ph.selectSpecies(f, 'H2')
        """

        if isinstance(species, types.StringTypes):
            k = self.speciesIndex(species)
            return f[k]
        elif species:
            fs = []
            k = 0
            for s in species:
                k = self.speciesIndex(s)
                fs.append(f[k])
            return asarray(fs)
        else:
            return asarray(f)

    def selectElements(self, f, elements):
        """Given an array *f* of floating-point element properties, return a
        those values corresponding to elements listed in *elements*. Returns an
        array if *elements* is a sequence, or a scalar if *elements* is a
        scalar.

        >>> f = ph.elementPotentials()
        >>> lam_o, lam_h = ph.selectElements(f, ['O', 'H'])
        >>> lam_h = ph.selectElements(f, 'H')
        """
        if isinstance(elements, types.StringTypes):
            m = self.elementIndex(elements)
            return f[m]
        if elements:
            fs = []
            k = 0
            for s in elements:
                k = self.elementIndex(s)
                fs.append(f[k])
            return asarray(fs)
        else:
            return asarray(f)
