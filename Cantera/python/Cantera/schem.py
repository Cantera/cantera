deprecated

"""
Surface chemistry.

Classes:

BulkSpecies   -- bulk species
SurfSpecies   -- surface species
SurfReaction  -- surface reactions
Interface     -- interfaces

"""


from Cantera import CanteraError
from Cantera import units
from Cantera.num import array
from Cantera import ctsurf, constants, SurfWriter
import math
import types

_rtypes = {'conc':0, 'cov':1, 'bar':2, 'site':3}
        

class SurfSpecies:
    """
    Surface species.

    A SurfSpecies object is required for every surface species.
    """
    
    def __init__(self, phase, symbol = '', elements = {}, size = 1):
        """
        Constructor arguments:

        phase -- the surface phase the species belongs to.
        symbol -- the symbol used in writing the reaction equation.
        elements -- a dictionary mapping element symbols to atom numbers.
        size -- number of surface sites occuppied by this species.

        Example:
        h_surf = SurfSpecies(phase = surf, symbol = 'H_s',
                             elements = {'H':1}, size = 1.0)
        
        """
        
        self.symbol = symbol
        self._nAtoms = {}
        for e in elements.keys():
            self._nAtoms[e] = elements[e]
        self.elements = elements.keys()        # list of element symbols
        self.size = size
        self.phase = phase
        self.index = len(self.phase.species)   
        self.phase.addSpecies(self)            # add it to the surface phase


    def nAtoms(self, el):
        """
        Number of atoms of element with symbol 'el'.
        """
        try:
            return self._nAtoms[el]
        except:
            return 0.0       # symbol not in the dictionary
        


class BulkSpecies:
    """
    Instances of class BulkSpecies represent the bulk-phase species in
    surface reactions. A 'bulk-phase' species is any species in a 3D
    phase, whether a gas, solid, or liquid. The species parameters are
    taken from the bulk-phase object.
    """

    def __init__(self, phase = None, symbol = ''):
        self.phase = phase
        self.symbol = symbol
        self.index = phase.speciesIndex(symbol)
        self._nAtoms = {}
        el = phase.elementNames()
        for m in range(len(el)):
            na = self.phase.nAtoms(self.index, m)
            if na <> 0:
                self._nAtoms[el[m]] = na
        self.weight = self.phase.molecularWeights()[self.index]
        self.elements = self._nAtoms.keys()  # list of element symbols

    # number of atoms
    def nAtoms(self, el):
        try:
            return self._nAtoms[el]
        except:
            return 0.0


class RateParam:
    """
    
    Reaction rate parameterizations appear in many different forms.
    Instances of this class are used to collect together a set of
    parameters for a particular parameterization type.

    The rate is computed using a generalized law of mass action, where
    the bulk (3D) species may be specified either by concentration or
    partial pressure, and empirical reaction orders may be assigned
    for each species. The surface species may be specified using
    either surface concentration, or coverage. Finally, the computed
    rate may be specified to be either a rate per unit area, or per
    site.
    
    Arguments:
    
       bulk -- one of ['conc', 'bar'].  Specifies how bulk-phase
       reactant abundances appear in the rate expression. Default:
       'conc'.  If 'bar' is specified, the ideal gas law is used to
       convert from partial pressure to concentration units when the
       reaction is added to the surface mechanism.

       surf -- one of ['conc', 'cov']. Specifies how surface reactant
       abundances appear in the rate expression. Default: 'conc'.

       rate -- one of ['per_area', 'per_site']. Specifies whether the
       rate expression computes the rate per unit area (default), or
       per site. These differ only by a factor of the total site
       density.

       Note that these values are used to convert the reaction
       parameterization into concentration form at the time the
       reaction is added to the surface mechanism. 
       
"""
    def __init__(self, bulk = 'conc', surf = 'conc', result = 'per_area'):
        self.bulk = bulk
        self.surf = surf
        self.result = result
        if not self.bulk in ['conc', 'bar']:
            raise CanteraError('unknown bulk phase quantity type: '+self.bulk)
        if not self.surf in ['conc', 'cov']:
            raise CanteraError('unknown surface phase quantity type: '
                               +self.surf)
        if not self.result in ['per_area', 'per_site']:
            raise CanteraError('unknown rate type: '+self.result)
        
    
    
                 
# A class for surface reactions.

class SurfReaction:
    """
    rate_type --- a triplet of strings specifying how the reaction rate is
    expressed. 
    """

    def __init__(self, phase = None, reactants = [], order = {}, products = [],
                 rate = None, stick = None,
                 type = None):

        self.reactants = reactants
        self.products = products
        self.order = order
        
        self.stick = stick
        self.rate = rate
        self.phase = phase
        if type:
            self.type = type
        else:
            self.type = RateParam(bulk = 'conc',
                                  surf = 'conc',
                                  result = 'per_area')

        self.phase.addReaction(self)


    
    def __repr__(self):
        rstr = ''
        for r in self.reactants:
            rstr += r.symbol + ' + '
        rstr = rstr[:-3] + ' => '
        for p in self.products:
            rstr += p.symbol + ' + '
        rstr = rstr[:-3]
        return rstr
        


class Interface:
    def __init__(self, phases = (None, None),
                 site_density = 0.0):

        self.rxns = []
        self.indexmap = {}

        # bulk phase 1 parameters
        self.p1 = phases[0]
        self.p1_nsp = 0
        self.p1id = 0
        self.p1sp = None
        if self.p1:
            self.p1id = self.p1.cthermo
            self.p1_nsp = self.p1.nSpecies()
            self.p1sp = self.p1.speciesNames()

        # bulk phase 2 parameters
        self.p2_nsp = 0
        self.p2id = 0
        self.p2 = phases[1]
        self.p2sp = None        
        if self.p2:
            self.p2id = self.p2.cthermo
            self.p2_nsp = self.p2.nSpecies()
            self.p2sp = self.p2.speciesNames()            

        self.s0 = site_density
        self.species = []
        self._spsymbols = []
        self._elements = {}
        self._freeze_species = 0
        self.rindex = []
        self.pindex = [] 

        # unit conversion factors to SI
        self._umol    = 1.0
        self._ulength = 1.0
        self._utime   = 1.0
        self._uconc2  = 1.0
        self._uconc3  = 1.0
        self._ue = 1.0

        # create instances of kernel classes SurfacePhase and
        # SurfKinetics
        self.__surf_id = ctsurf.surf_new(self.s0)
        self.__surfkin_id = ctsurf.surfkin_new(self.__surf_id,
                                               self.p1id, self.p2id)


        self.rxndata = []

        
    def __del__(self):
        ctsurf.surf_delete(self.__surf_id)
        ctsurf.surfkin_delete(self.__surfkin_id)
        ctsurf.surf1d_delete(self.__surf1d_id)


    def kin_id(self):
        return self.__surfkin_id

    def surf_id(self):
        return self.__surf1d_id

##     def parse_equation(self, eqn):
##         toks = string.split(eqn)
##         digits = ['2', '3', '4', '5', '6', '7', '8', '9']
##         for t in toks:
##             if t in digits:
##                 nu = int(t)
##             elif t == '+'
                
##             if t in self.p1sp:
##                 k = self.p1sp.index(t)
##             elif t in self.p2sp:
##                 k = self.p1nsp + self.p2sp.index(t)
##             elif t in self._spsymbols:
##                 k = self.p1nsp + self.p2nsp + self._spsymbols.index(t)
##             else:
##                 raise CanteraError('unknown species '+t)
        

    def addDoc(self, key, value):
        ctsurf.surf_doc(self.__surf_id, key, value)

    def save(self, filename, idtag, comment):
        ctsurf.surfkin_save(self.__surfkin_id, filename, idtag, comment)
        
    def setUnits(self, 
                 length = '', moles = '', time = '', E = ''):
        
        """
        Set the unit system for surface reactions. On each call,
        only the specified unit conversion factors are set; the rest
        are left unchanged.

        length -- units for lengths ('m', 'cm', 'mm')
        moles  -- units for amount ('mol', 'kmol', 'molecule')
        time   -- units for time ('s')
        E      -- units for activation energy ('kcal_per_mol',
                  'cal_per_mol', 'eV', 'K')
        """
        if moles: self._umol    = units.mole(moles)
        if length: self._ulength = units.length(length)
        if time: self._utime   = 1.0
        if E: self._ue = units.actEnergy(E)
        self._uconc2 = self._umol/(self._ulength * self._ulength)
        self._uconc3 = self._uconc2/self._ulength

        
    def addSpecies(self, s):
        """Add a species."""
        if self._freeze_species > 0:
            raise CanteraError("Species must be added to a phase "+
                               "before reactions are added")
        self.species.append(s)
        self._spsymbols.append(s.symbol)
        self.indexmap[s.symbol] = s.index
        for e in s.elements:
            self._elements[e] = 1
        self.elements = self._elements.keys()
        ctsurf.surf_addspecies(self.__surf_id, s.symbol, s.size)

    def nSpecies(self):
        return ctsurf.surf_nspecies(self.__surf_id)

    def speciesIndex(self, sym):
        return self.indexmap[sym]

    def speciesNames(self):
        return self._spsymbols
    
    def setSiteDensity(self, s0):
        self.s0 = s0
        print 'setting site density to ',s0
        ctsurf.surf_setsitedensity(self.__surf_id, self.s0)

    def show(self):
        cv = self.coverages()
        c = map(None,self._spsymbols,cv)        
        for s in c:
            print '  %16s     %10.4g ' % (s[0],s[1])
            
    def coverage(self, sp):
        clist = ctsurf.surf_getcoverages(self.__surf_id)
        c = []
        for s in sp:
            c.append(clist[self.indexmap[s]])
        return c
    
    def coverages(self):
        return ctsurf.surf_getcoverages(self.__surf_id)

    def setCoverages(self, cov):
        if type(cov) == types.DictType:
            cv = [0.0]*len(self.species)
            for c in cov.keys():
                cv[self.indexmap[c]] = cov[c]
        else:
            cv = cov
        ctsurf.surf_setcoverages(self.__surf_id, array(cv, 'd'))

    def ratesOfProgress(self):
        """Rates of progress for all surface reactions [kmol/m^2-s]."""
        return ctsurf.surfkin_getratesofprogress(self.__surf_id)

    def netProductionRates(self):
        """Net production rates for all bulk and surface species.

        The results are returned as a tuple of 3 arrays for the species in
        bulk phase 1, bulk phase 2, and the surface phase, respectively.
        """
        n1 = self.p1_nsp
        n2 = self.p2_nsp
        ns = self.nSpecies()
        ntot = n1 + n2 + ns   
        sdot = ctsurf.surfkin_getsdot(self.__surf_id, ntot)
        return (sdot[:n1], sdot[n1:n1+n2], sdot[n1+n2:])

    def integrate(self, dt):
        """Integrate the surface site-balance equations for time dt
        with fixed bulk compositions."""
        ctsurf.surfkin_integrate(self.__surf_id, dt)

    def _index(self, sp):
        """Location of species object 'sp' in the solution array.""" 
        if sp.phase == self.p1:
            return sp.index
        elif sp.phase == self.p2:
            return sp.index + self.p1_nsp
        elif sp.phase == self:
            return sp.index + self.p1_nsp + self.p2_nsp

    def _uconc(self, n):
        """Return the concentration conversion factor to SI for the
        species with index n."""
        if n < self.p1_nsp + self.p2_nsp:
            return self._uconc3
        else:
            return self._uconc2

    def isBulkSpecies(self, n):
        """Return 1 if species n is a bulk species, and 0 if it is a
        surface species."""
        if n < self.p1_nsp + self.p2_nsp:
            return 1
        else:
            return 0
        
    def _bulkReactants(self, reactants):
        """Given a list of reactant objects, return a list of those that
        are bulk species."""
        rb = []
        for r in reactants:
            if hasattr(r,'weight'):
                rb.append(r)
        return rb

    def _surfReactants(self, reactants):
        """Given a list of reactant objects, return a list of those that
        are surface species."""        
        rs = []
        for r in reactants:
            if not hasattr(r,'weight'):
                rs.append(r)
        return rs


    def write_C(self, mech):

        writer = SurfWriter.SurfWriter(self, mech)
        for d in self.rxndata:
            rindex, rstoich, rorder, pindex, pstoich, rate = d
            writer.write_ROP(rindex, rstoich, rorder,
                             pindex, pstoich, rate) 
            writer.write_sdot(rindex, rstoich, pindex, pstoich)
            writer.write_update_rate(rate)
            writer.nrxns += 1
        writer.write()

        
    def addReaction(self, r):
        """Add a reaction to the surface mechanism."""
        self._freeze_species = 1
            
        self.rxns.append(r)

            
        # check balance
        for e in self.elements:
            n = 0
            for rr in r.reactants:
                n += rr.nAtoms(e)
            for pp in r.products:
                n -= pp.nAtoms(e)
            if n <> 0:
                raise CanteraError('Reaction does not balance. \nElement '+
                                   e+' out of balance by '+`n`+' atoms.')

        rmap = {}
        rindex2object = {}
        for rr in r.reactants:
            if rmap.has_key(rr):
                rmap[rr][1] += 1
                rmap[rr][2] += 1                
            else:
                rmap[rr] = [self._index(rr),1,1]
                rindex2object[self._index(rr)] = rr

        pmap = {}
        for pp in r.products:
            if pmap.has_key(pp):
                pmap[pp][1] += 1
                pmap[pp][2] += 1                
            else:
                pmap[pp] = [self._index(pp),1,1]            

        rinfo = array(rmap.values(),'i')
        pinfo = array(pmap.values(),'i')

        rindex = rinfo[:,0]
        for rr in r.order.keys():
            loc = self._index(rr)
            for n in range(len(rindex)):
                if loc == rindex[n]:
                    rinfo[n,2] = r.order[rr]
    
        rstoich = rinfo[:,1]
        rorder  = rinfo[:,2]
        pindex  = pinfo[:,0]
        pstoich = pinfo[:,1]
    
        funit = 1.0
        fb = 0.0


        # sticking probabilities are dimensionless, and are limited to
        # reactions with one bulk reactant. The reaction rate is
        # computed as the sticking probability multiplied by the
        # incident flux, times the coverages of all surface species
        # reactants.
        
        if r.stick:
            rb = self._bulkReactants(r.reactants)            
            if len(rb) <> 1:
                raise CanteraError(
                    'A sticking probability can only be specified if there'+
                    ' is exactly one bulk-phase reactant')
            wt = rb[0].weight
            
            # mean speed / sqrt(T)
            cfactor = math.sqrt(8.0*constants.GasConstant/(constants.Pi * wt))
            
            for n in range(len(rindex)):
                if not self.isBulkSpecies(rindex[n]):
                    sz = rindex2object[rindex[n]].size
                    funit /= math.pow(self.s0/sz, rorder[n])
            preExp = funit*r.stick[0]*cfactor/4.0
            b = r.stick[1] + 0.5
            e = r.stick[2]
            r.rate = [preExp, b, e]


        # Reaction rates are specified using several different
        # conventions.  Here the various possibilities are converted
        # into rates per unit area, with reactant amounts specified in
        # concentrations.
        
        else:
            for n in range(len(rindex)):
                if self.isBulkSpecies(rindex[n]):
                    if r.type.bulk == 'conc':
                        funit /= (math.pow(self._uconc3,rorder[n]))
                    elif r.type.bulk == 'bar':
                        funit *= math.pow(1.e-5
                                          * constants.GasConstant,rorder[n])
                        fb += rorder[n]
                else:
                    if r.type.surf == 'conc':
                        funit /= (math.pow(self._uconc2,rorder[n]))
                    elif r.type.surf == 'cov':
                        sz = rindex2object[rindex[n]].size
                        funit /= math.pow(self.s0/sz, rorder[n])
                    
            if r.type.result == 'per_area':
                funit *= self._uconc2
            elif r.type.result == 'per_site':
                print 'multiplying A by ',self.s0
                funit *= self.s0

            # modify the pre-exponential term, and the temperature
            # exponent
            r.rate[0] *= funit
            print 'A = ',r.rate[0] 
            r.rate[1] += fb

            
        rate = array(r.rate,'d')

        # convert the activation energy to K
        rate[2] *= self._ue

        self.rxndata.append((rindex, rstoich, rorder, pindex, pstoich, rate))
        
        ctsurf.surfkin_addreaction(self.__surf_id,
                                   len(rindex),
                                   array(rindex,'i'),
                                   array(rstoich,'i'),
                                   array(rorder,'i'),
                                   len(pindex),
                                   array(pindex,'i'),
                                   array(pstoich,'i'),
                                   len(r.rate), rate)
               
