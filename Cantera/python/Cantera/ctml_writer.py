
#
# Cantera input file processor
#

from Cantera import GasConstant
from Cantera.XML import XML_Node
import types, math

SPECIES = 10
SPECIES_SET = 20
COLLECTION = 30
THERMO = 40

ALL = 10
SKIP_UNDECLARED_ELEMENTS = 20
SKIP_UNDECLARED_SPECIES = 30
STOP = 0

# default units
_ulen = 'm' 
_umol = 'kmol'
_umass = 'kg'
_utime = 's'
_ue = 'J/kmol'

# default std state pressure
_pref = 1.0e5    # 1 bar

_name = 'noname'

_species = []
_phases = []
_reactions = []
_atw = {}
_mw = {}

def isnum(a):
    if type(a) == types.IntType or type(a) == types.FloatType:
        return 1
    else:
        return 0
    
def dataset(nm):
    global _name
    _name = nm

def standard_pressure(p0):
    """Set the default standard-state pressure."""
    global _pref
    _pref = p0
    
def get_atomic_wts():
    global _atw
    edb = XML_Node('edb', src = 'elements.xml')
    edata = edb.child('ctml/elementData')
    e = edata.children()
    for el in e:
        if el['name'] <> 'dummy':
            _atw[el['name']] = el['atomicWt']
            if el['atomicWt'] == '':
                print 'no atomic weight for ',el['name']
    
def units(length = '', quantity = '', mass = '', time = '', act_energy = ''):
    global _ulen, _umol, _ue
    if length: _ulen = length
    if quantity: _umol = quantity
    if act_energy: _ue = act_energy
    if time: _utime = time
    if mass: _umass = mass

def ufmt(base, n):
    if n == 0: return ''
    if n == 1: return '-'+base
    if n == -1: return '/'+base
    if n > 0: return '-'+base+`n`
    if n < 0: return '/'+base+`-n`    
    
def write():
    x = XML_Node("ctml")
    for ph in _phases:
        ph.build(x)
    s = species_set(name = _name, species = _species)
    s.build(x)
    
    r = x.addChild('reactionData')
    r['id'] = 'reaction_data'
    for rx in _reactions:
        rx.build(r)
    
    if _name <> 'noname':
        x.write(_name+'.xml')
    else:
        print x
            
def addFloat(x, nm, val, fmt=''):
    u = ''
    s = ''
    if isnum(val):
        fval = float(val)
        if fmt:
            s = fmt % fval
        else:
            s = `fval`
        x.addChild(nm, s)
    else:
        v = val[0]
        u = val[1]
        if fmt:
            s = fmt % v
        else:
            s = `v`
        xc = x.addChild(nm, s)
        xc['units'] = u
    
def getAtomicComp(atoms):
    if type(atoms) == types.DictType: return atoms
    a = atoms.replace(',',' ')
    toks = a.split()
    d = {}
    for t in toks:
        b = t.split(':')
        d[b[0]] = int(b[1])
    return d

def getEfficiencies(e):
    if type(e) == types.DictType: return e
    a = atoms.replace(',',' ')
    toks = a.split()
    d = {}
    for t in toks:
        b = t.split(':')
        d[b[0]] = int(b[1])
    return d

def getReactionSpecies(s):
    toks = s.replace(' + ',' ').split()
    d = {}
    n = 1
    for t in toks:
        if t > '0' and t < '9':
            n = int(t)
        else:
            if d.has_key(t):
                d[t] += n
            else:
                d[t] = n
            n = 1
    return d

class writer:
    def write_ctml(self, file = ''):
        x = XML_Node("ctml")
        self.build(x)
        if file:
            x.write(file)
        else:
            print x
    
class collection(writer):
    def __init__(self, s):
        self._s = s
        self.type = COLLECTION        
    def build(self, p):
        for s in self._s:
            s.build(p)
            
class species_set(writer):
    def __init__(self, name = '', species = []):
        self._s = species
        self._name = name
        self.type = SPECIES_SET
        
    def build(self, p):
        p.addComment('     species definitions     ')
        sd = p.addChild("speciesData")
        sd.addAttrib("id","species_data")
        for s in self._s:
            if s.type == SPECIES:
                s.build(sd)
            else:
                raise 'wrong object type in species_set: '+s.__class__

            
class species(writer):
    """A species."""
    
    def __init__(self,
                 name = 'missing name!',
                 atoms = '',
                 comment = '',
                 thermo = None,
                 transport = None,
                 charge = -999):
        self._name = name
        self._atoms = getAtomicComp(atoms)
        mw = 0.0
        for a in self._atoms.keys():
            mw += self._atoms[a]*float(_atw[a])
        self._mw = mw
        global _mw
        _mw[name] = mw
        self._comment = comment
        self._thermo = thermo
        self._transport = transport
        chrg = 0
        self._charge = charge
        if self._atoms.has_key('E'):
            chrg = -self._atoms['E']
            if self._charge <> -999:
                if self._charge <> chrg:
                    raise 'specified charge inconsistent with number of electrons'
            else:
                self._charge = chrg
        self.type = SPECIES
        
        global _species
        _species.append(self)
        

    def build(self, p):
        #phname = ''
        #for ph in _phases:
        #    if ph.has_species(self._name):
        #        phname = ph._name
        hdr = '    species '+self._name+'    '
        p.addComment(hdr)        
        s = p.addChild("species")
        s.addAttrib("name",self._name)
        a = ''
        for e in self._atoms.keys():
            a += e+':'+`self._atoms[e]`+' '
        s.addChild("atomArray",a)
        if self._comment:
            s.addChild("note",self._comment)
        if self._charge <> -999:
            s.addChild("charge",self._charge)
        if self._thermo:
            t = s.addChild("thermo")
            if type(self._thermo) == types.InstanceType:
                self._thermo.build(t)
            else:
                nt = len(self._thermo)
                for n in range(nt):
                    self._thermo[n].build(t)
        if self._transport:
            t = s.addChild("transport")
            if type(self._transport) == types.InstanceType:
                self._transport.build(t)
            else:
                nt = len(self._transport)
                for n in range(nt):
                    self._transport[n].build(t)                    

class thermo(writer):
    """Base class for species standard-state thermodynamic properties."""
    def _build(self, p):
        return p.addChild("thermo")

        
class NASA(thermo):
    """NASA polynomial parameterization."""
    
    def __init__(self, range = (0.0, 0.0), 
                 coeffs = [], p0 = -1.0):
        self._t = range
        self._pref = p0
        if len(coeffs) <> 7:
            raise 'NASA coefficient list must have length = 7'
        self._coeffs = coeffs

        
    def build(self, t):
        n = t.addChild("NASA")
        n['Tmin'] = `self._t[0]`
        #n['Tmid'] = `self._t[1]`        
        n['Tmax'] = `self._t[1]`
        if self._pref <= 0.0:
            n['P0'] = `_pref`
        else:
            n['P0'] = `self._pref`
        str = ''
        for i in  range(7):
            str += '%17.9E, ' % self._coeffs[i]
            if i > 0 and 3*((i+1)/3) == i: str += '\n'
        str = str[:-2]
        u = n.addChild("floatArray", str)
        u["size"] = "7"
        u["name"] = "coeffs"
    
                 
class const_cp(thermo):
    """Constant specific heat."""
    
    def __init__(self, tmax = -1.0, tmin = -1.0,
                 t0 = 298.15, cp0 = 0.0, h0 = 0.0, s0 = 0.0):
        self._t = [tmin, tmax]
        self._c = [t0, h0, s0, cp0]
        
    def build(self, t):
        #t = self._build(p)
        c = t.addChild('const_cp')
        if self._t[0] >= 0.0: c['Tmin'] = self._t[0]
        if self._t[1] >= 0.0: c['Tmax'] = self._t[1]
        addFloat(c,'t0',self._c[0])
        addFloat(c,'h0',self._c[1])
        addFloat(c,'s0',self._c[2])        
        addFloat(c,'cp0',self._c[3])


class gas_transport:
    """Transport coefficients for ideal gas transport model."""

    def __init__(self, geom = 'nonlin',
                 diam = 0.0, well_depth = 0.0, dipole = 0.0,
                 polar = 0.0, rot_relax = 0.0):
        self._geom = geom
        self._diam = diam
        self._well_depth = well_depth
        self._dipole = dipole
        self._polar = polar
        self._rot_relax = rot_relax

    def build(self, t):
        #t = s.addChild("transport")
        t['model'] = 'gas_transport'
        #        t.addChild("geometry", self._geom)
        tg = t.addChild('string',self._geom)
        tg['title'] = 'geometry'
        addFloat(t, "LJ_welldepth", (self._well_depth, 'K'), '%8.3f')
        addFloat(t, "LJ_diameter", (self._diam, 'A'),'%8.3f')
        addFloat(t, "dipoleMoment", (self._dipole, 'Debye'),'%8.3f')        
        addFloat(t, "polarizability", (self._polar, 'A3'),'%8.3f')
        addFloat(t, "rotRelax", self._rot_relax,'%8.3f')        
        
class Arrhenius(writer):
    def __init__(self,
                 A = 0.0,
                 n = 0.0,
                 E = 0.0):
        self._c = [A, n, E]
        
    def build(self, p, units = '', gas_species = [], name = ''):
        a = p.addChild('Arrhenius')
        if name: a['name'] = name
        if isnum(self._c[0]):
            addFloat(a,'A', (self._c[0], units), fmt = '%14.6E')
        else:
            addFloat(a,'A',self._c[0], fmt = '%14.6E')
        a.addChild('b',`self._c[1]`)
        if isnum(self._c[2]):
            addFloat(a,'E',(self._c[2],_ue), fmt = '%f')
        else:
            addFloat(a,'E',self._c[2], fmt = '%f')
                     

class stick(writer):
    def __init__(self,
                 A = 0.0,
                 n = 0.0,
                 E = 0.0):
        self._c = [A, n, E]
        self._sp = species
        
    def build(self, p, units = '', gas_species = [], name = ''):
        a = p.addChild('Stick')
        ngas = len(gas_species)
        if ngas <> 1:
                raise 'sticking probabilities can only be used for reactions with one gas-phase reactant'
        a['species'] = gas_species[0]
        if name: a['name'] = name
        addFloat(a,'A',self._c[0], fmt = '%14.6E')
        a.addChild('b',`self._c[1]`)
        if isnum(self._c[2]):
            addFloat(a,'E',(self._c[2],_ue), fmt = '%f')
        else:
            addFloat(a,'E',self._c[2], fmt = '%f')

        
class reaction(writer):
    def __init__(self,
                 equation = '',
                 kf = None,
                 id = '',                 
                 ):

        self._id = id
        self._e = equation
        global _reactions
        self._num = len(_reactions)+1
        r = ''
        p = ''
        for e in ['<=>', '=>', '=']:
            if self._e.find(e) >= 0:
                r, p = self._e.split(e)
                if e in ['<=>','=']: self.rev = 1
                else: self.rev = 0
                break
        self._r = getReactionSpecies(r)
        self._p = getReactionSpecies(p)
        self._kf = kf
        self._igspecies = []
        self._type = ''
        _reactions.append(self)

        
    def build(self, p):
        if self._id:
            id = self._id
        else:
            if self._num < 10:
                nstr = '000'+`self._num`
            elif self._num < 100:
                nstr = '00'+`self._num`
            elif self._num < 1000:
                nstr = '0'+`self._num`                
            else:
                nstr = `self._num`
            id = nstr
        p.addComment("   reaction "+id+"    ")                
        r = p.addChild('reaction')
        r['id'] = id
        if self.rev:
            r['reversible'] = 'yes'
        else:
            r['reversible'] = 'no'

        ee = self._e.replace('<','[')
        ee = ee.replace('>',']')
        r.addChild('equation',ee)
        
        mdim = 0
        ldim = 0
        str = ''

        for s in self._r.keys():
            ns = self._r[s]
            nm = -999
            nl = -999

            str += s+':'+`ns`+' '
            
            for ph in _phases:
                if ph.has_species(s):
                    nm, nl = ph.conc_dim()
                    if ph.is_ideal_gas():
                        self._igspecies.append(s)
                    break
            if nm < 0:
                print self._r
                raise 'undeclared species '+s
            else:
                mdim += nm*ns
                ldim += nl*ns

        # adjust the moles and length powers based on the dimensions of
        # the rate of progress (moles/length^2 or moles/length^3)
        if self._type == 'surface':
            mdim += -1
            ldim += 2
        else:
            mdim += -1
            ldim += 3


        # add the reaction type as an attribute if it has been specified.
        if self._type:
            r['type'] = self._type

        # The default rate coefficient type is Arrhenius. If the rate
        # coefficient has been specified as a sequence of three
        # numbers, then create a new Arrhenius instance for it;
        # otherwise, just use the supplied instance.
        nm = ''
        kfnode = r.addChild('rateCoeff')
        if self._type == '':
            self._kf = [self._kf]
        elif self._type == 'surface':
            self._kf = [self._kf]            
        elif self._type == 'threeBody':
            self._kf = [self._kf]
            mdim += 1
            ldim -= 3
            
        for kf in self._kf:
            
            # compute the pre-exponential units string, and if it begins with a
            # dash, remove it.
            ku = ufmt(_ulen,-ldim) + ufmt(_umol,-mdim) + ufmt('s',-1)
            if ku[0] == '-': ku = ku[1:]

            if type(kf) == types.InstanceType:
                k = kf
            else:
                k = Arrhenius(A = kf[0], n = kf[1], E = kf[2])
            k.build(kfnode, units = ku, gas_species = self._igspecies, name = nm)

            # set values for low-pressure rate coeff if falloff rxn
            mdim += 1
            ldim -= 3
            nm = 'k0'

        str = str[:-1]
        r.addChild('reactants',str)
        str = ''
        for s in self._p.keys():
            ns = self._p[s]
            str += s+':'+`ns`+' '
        str = str[:-1]
        r.addChild('products',str)
        return r

#-------------------

class three_body_reaction(reaction):
    def __init__(self,
                 equation = '',
                 kf = None,
                 efficiencies = '',
                 id = '',                 
                 ):

        reaction.__init__(self, equation, kf, id)
        self._type = 'threeBody'
        self._effm = 1.0
        self._eff = efficiencies

        # clean up reactant and product lists
        for r in self._r.keys():
            if r == 'M' or r == 'm':
                del self._r[r]
        for p in self._p.keys():
            if p == 'M' or p == 'm':
                del self._p[p]
                

        
    def build(self, p):
        r = reaction.build(self, p)
        kfnode = r.child('rateCoeff')
        
        if self._eff:
            eff = kfnode.addChild('efficiencies',self._eff)
            eff['default'] = `self._effm`
        

#---------------


class falloff_reaction(reaction):
    def __init__(self,
                 equation = '',
                 kf0 = None,
                 kf = None, 
                 efficiencies = '',
                 falloff = None,
                 id = ''
                 ):

        kf2 = (kf, kf0)        
        reaction.__init__(self, equation, kf2, id)
        self._type = 'falloff'
        # use a Lindemann falloff function by default
        self._falloff = falloff
        if self._falloff == None:
            self._falloff = Lindemann()
        
        self._effm = 1.0
        self._eff = efficiencies        

        # clean up reactant and product lists
        del self._r['(+']
        del self._p['(+']
        if self._r.has_key('M)'):
            del self._r['M)']
            del self._p['M)']
        if self._r.has_key('m)'):
            del self._r['m)']
            del self._p['m)']
        else:
            for r in self._r.keys():
                if r[-1] == ')' and r.find('(') < 0:
                    if self._eff:
                        raise '(+ '+mspecies+') and '+self._eff+' cannot both be specified'
                    self._eff = r[-1]+':1.0'
                    self._effm = 0.0
                    
                    del self._r[r]
                    del self._p[r]

        
    def build(self, p):
        r = reaction.build(self, p)

        kfnode = r.child('rateCoeff')
        
        if self._eff and self._effm >= 0.0:
            eff = kfnode.addChild('efficiencies',self._eff)
            eff['default'] = `self._effm`
        
        if self._falloff:
            self._falloff.build(kfnode)


class surface_reaction(reaction):

    def __init__(self,
                 equation = '',
                 kf = None,
                 stick = None,
                 id = ''):
        if stick:
            reaction.__init__(self, equation, stick, id)
        else:
            reaction.__init__(self, equation, kf, id)
        self._type = 'surface'

#--------------            


class state:
    def __init__(self,
                 temperature = None,
                 pressure = None,
                 mole_fractions = None,
                 mass_fractions = None,
                 density = None,
                 coverages = None):
        self._t = temperature
        self._p = pressure
        self._rho = density
        self._x = mole_fractions
        self._y = mass_fractions
        self._c = coverages
        
    def build(self, ph):
        st = ph.addChild('state')
        if self._t: addFloat(st, 'temperature', self._t)
        if self._p: addFloat(st, 'pressure', self._p)
        if self._rho: addFloat(st, 'density', self._rho)
        if self._x: st.addChild('moleFractions', self._x)
        if self._y: st.addChild('massFractions', self._y)
        if self._c: st.addChild('coverages', self._c)
    

class phase(writer):
    """Base class for phases of matter."""
    
    def __init__(self,
                 name = '',
                 dim = 3,
                 elements = '',
                 species = '',
                 reactions = 'none',
                 initial_state = None):
        
        self._name = name
        self._dim = dim
        self._el = elements
        self._sp = []
        self._rx = []

        #--------------------------------
        #        process species
        #--------------------------------
        
        # if a single string is entered, make it a list
        if type(species) == types.StringType:
            self._species = [species]
        else:
            self._species = species

        self._skip = 0

        # dictionary of species names 
        self._spmap = {}

        # for each species string, check whether or not the species
        # are imported or defined locally. If imported, the string
        # contains a colon (:) 
        for sp in self._species:
            if sp.find(':') > 0:
                datasrc, spnames = sp.split(':')
                self._sp.append((datasrc+'.xml', spnames))
            else:
                spnames = sp
                self._sp.append(('', spnames))                

            # strip the commas, and make the list of species names
            sptoks = spnames.replace(',',' ').split()            
            for s in sptoks:
                if self._spmap.has_key(s):
                    raise 'multiply-defined species '+s+' in phase '+self._name
                self._spmap[s] = self._dim
            
        self._rxns = reactions

        # check that species have been declared
        if len(self._spmap) == 0:
            raise 'No species declared for phase '+self._name
        
        # and that only one species is declared if it is a pure phase
        if self.is_pure() and len(self._spmap) > 1:
            raise 'Pure phases may only declare one species, but phase '+self._name+' declares '+`len(self._spmap)`+'.'

        self._initial = initial_state

        # add this phase to the global phase list
        global _phases
        _phases.append(self)

    def is_ideal_gas(self):
        return 0
    
    def is_pure(self):
        return 0
    
    def has_species(self, s):
        """Return 1 is a species with name 's' belongs to the phase,
        or 0 otherwise."""
        if self._spmap.has_key(s): return 1
        return 0

    def conc_dim(self):
        """Concentration dimensions. Used in computing the units for reaction
        rate coefficients."""
        return (1, -self._dim) 

    
    def buildrxns(self, p):
        
        if type(self._rxns) == types.StringType:
            self._rxns = [self._rxns]
            
        # for each reaction string, check whether or not the reactions
        # are imported or defined locally. If imported, the string
        # contains a colon (:) 
        for r in self._rxns:
            
            if r.find(':') > 0:
                datasrc, rnum = r.split(':')
                self._rx.append((datasrc+'.xml', rnum))
            else:
                rnum = r
                self._rx.append(('', rnum))                

        for r in self._rx:
            datasrc = r[0]
            ra = p.addChild('reactionArray')
            ra['datasrc'] = datasrc+'#reaction_data'
            rtoks = r[1].split()
            if rtoks[0] <> 'all':
                i = ra.addChild('include')
                #i['prefix'] = 'reaction_'
                i['min'] = rtoks[0]
                if len(rtoks) > 2 and (rtoks[1] == 'to' or rtoks[1] == '-'):
                    i['max'] = rtoks[2]
                else:
                    i['max'] = rtoks[0]                
                    
            
    def build(self, p):
        p.addComment('    phase '+self._name+'     ')
        ph = p.addChild('phase')
        ph['id'] = self._name
        ph['dim'] = `self._dim`
        e = ph.addChild('elementArray',self._el)
        e['datasrc'] = 'elements.xml'
        for s in self._sp:
            datasrc, names = s
            sa = ph.addChild('speciesArray',names)
            sa['datasrc'] = datasrc+'#species_data'
            
        if self._skip:
            sk = sa.addChild('skip')
            sk['element'] = 'undeclared'
            
        if self._rxns <> 'none':
            self.buildrxns(ph)
            
        #self._eos.build(ph)
        if self._initial:
            self._initial.build(ph)
        return ph


class ideal_gas(phase):
    """An ideal gas mixture."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 reactions = 'none',
                 kinetics = 'GasKinetics',
                 transport = 'None',
                 initial_state = None):
        
        phase.__init__(self, name, 3, elements, species, reactions,
                       initial_state)
        self._pure = 0
        self._kin = kinetics
        self._tr = transport

        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'IdealGas'
        k = ph.addChild("kinetics")
        k['model'] = self._kin
        t = ph.addChild('transport')
        t['model'] = self._tr

    def is_ideal_gas(self):
        return 1
    
class pure_solid(phase):
    """A pure solid."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None):
        
        phase.__init__(self, name, 3, elements, species, 'none',
                       initial_state)
        self._dens = density
        self._pure = 1
        if self._dens < 0.0:
            raise 'density must be specified.'        
        self._pure = 0
        self._tr = transport

    def conc_dim(self):
        return (0,0)
        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'SolidCompound'
        addFloat(e, 'density', self._dens)
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'        


class ideal_interface(phase):
    """An ideal interface."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 reactions = 'none',
                 site_density = 0.0,
                 phases = [],
                 kinetics = 'Interface',
                 transport = 'None',
                 initial_state = None):

        self._type = 'surface'
        phase.__init__(self, name, 2, elements, species, reactions,
                       initial_state)
        self._pure = 0
        self._kin = kinetics
        self._tr = transport
        self._phases = phases
        self._sitedens = site_density
        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'Surface'
        addFloat(e, 'site_density',self._sitedens)
        k = ph.addChild("kinetics")
        k['model'] = self._kin
        t = ph.addChild('transport')
        t['model'] = self._tr
        p = ph.addChild('phaseArray',self._phases)


    def conc_dim(self):
        return (1, -2)

        
#------------------ equations of state --------------------------

class eos(writer):
    def is_pure(self):
        return self._pure
    
class incompressible_eos(eos):
    def __init__(self, density = -1.0):
        self._dens = density
        self._pure = 0
        if self._dens < 0.0:
            raise 'density must be specified.'
        
    def build(self, p):
        e = p.addChild("thermo")
        e['model'] = 'Incompressible'
        addFloat(e, 'density', self._dens)

    def conc_dim(self):
        return (1, -3)
    
class solid_compound_eos(eos):
    def __init__(self, density = -1.0):
        self._dens = density
        self._pure = 1
        if self._dens < 0.0:
            raise 'density must be specified.'
        
    def build(self, p):
        e = p.addChild("thermo")
        e['model'] = 'SolidCompound'
        addFloat(e, 'density', self._dens)
        if len(self.parent._spmap) > 1:
            raise 'A solid compound can only have one species.'
        
    def conc_dim(self):
        return (0, 0)


class ideal_gas_eos(eos):
    def __init__(self, kinetics = 'GasKinetics',
                 transport = 'none'):
        self._pure = 0
        self._kin = kinetics
        self._tr = transport
        global _idealgas_class
        _idealgas_class = self.__class__
        
    def build(self, p):
        e = p.addChild("thermo")
        e['model'] = 'IdealGas'
        k = p.addChild("kinetics")
        k['model'] = self._kin
        t = p.addChild('transport')
        t['model'] = self._tr
        
    def conc_dim(self):
        return (1, -3)


class surface(eos):
    def __init__(self, site_density = 0.0):
        self._pure = 0
        self._s0 = site_density
    def build(self, p):
        e = p.addChild("thermo")
        e['model'] = 'Surface'
        addFloat(e, 'site_density', self._s0, '%14.6E')
        
    def conc_dim(self):
        return (1, -2)    
        
#-------------------------------------------------------------------

# falloff parameterizations

class Troe:
    
    def __init__(self, A = 0.0, T3 = 0.0, T1 = 0.0, T2 = -999.9):
        if T2 <> -999.9:
            self._c = (A, T3, T1, T2)
        else:
            self._c = (A, T3, T1)
            
    def build(self, p):
        s = ''
        for num in self._c:
            s += '%g ' % num
        f = p.addChild('falloff', s)
        f['type'] = 'Troe'


class SRI:
    def __init__(self, A = 0.0, B = 0.0, C = 0.0, D = -999.9, E=-999.9):
        if D <> -999.9 and E <> -999.9:
            self._c = (A, B, C, D, E)
        else:
            self._c = (A, B, C)
        
    def build(self, p):
        s = ''
        for num in self._c:
            s += '%g ' % num
        f = p.addChild('falloff', s)
        f['type'] = 'SRI'

class Lindemann:
    def __init__(self):
        pass
    def build(self, p):
        f = p.addChild('falloff')
        f['type'] = 'Lindemann'

#--------------------------------------------------------------------

## class gas_transport:
##     def __init__(self, geom = 'nonlin',
##                  welldepth = 0.0,
##                  diam = 0.0,
##                  dipole = 0.0,
##                  polar = 0.0,
##                  rot_relax = 0.0):
##         self._sp = species
##         self._geom = geom
##         self._params = (welldepth, diam, dipole, polar, rotrelax)

##         #global _trdata
##         #_trdata[species] = self

##     def build(self, s):
##         tr = s.addChild('transport')
##         g = tr.addChild('string','linear')
##         g['title'] = 'geometry'
##         tr.addChild('LJ_welldepth',`self._params[0]`)
##         tr.addChild('LJ_diameter',`self._params[1]`)
##         tr.addChild('dipoleMoment',`self._params[2]`)
##         tr.addChild('polarizability',`self._params[3]`)
##         tr.addChild('rotRelax',`self._params[4]`)        


get_atomic_wts()
