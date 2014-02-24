#!/usr/bin/env python

##
# @file ctml_writer.py
#
# Cantera .cti input file processor
# @defgroup pygroup Cantera Python Interface
#
# The functions and classes in this module process Cantera .cti input
# files and produce CTML files. It can be imported as a module, or used
# as a script.
#
# script usage:
#
# python ctml_writer.py infile.cti
#
# This will produce CTML file 'infile.xml'

from __future__ import print_function

class CTI_Error(Exception):
    """Exception raised if an error is encountered while
    parsing the input file.
    @ingroup pygroup"""
    def __init__(self, msg):
        print('\n\n***** Error parsing input file *****\n\n')
        print(msg)
        print()




indent = ['',
          ' ',
          '  ',
          '   ',
          '    ',
          '     ',
          '      ',
          '       ',
          '        ',
          '          ',
          '           ',
          '            ',
          '             ',
          '              ',
          '               ',
          '                ']

#-----------------------------------------------------

class XMLnode(object):
    """This is a minimal class to allow easy creation of an XML tree
    from Python. It can write XML, but cannot read it."""

    __slots__ = ('_name', '_value', '_attribs', '_children', '_childmap')

    def __init__(self, name="--", value = ""):

        """Create a new node. Usually this only needs to be explicitly
        called to create the root element. Method addChild calls this
        constructor to create the new child node."""

        self._name = name

        # convert 'value' to a string if it is not already, and
        # strip leading whitespace
        if not isinstance(value, str):
            self._value = repr(value).lstrip()
        else:
            self._value = value.lstrip()

        self._attribs = {}    # dictionary of attributes
        self._children = []   # list of child nodes
        self._childmap = {}   # dictionary of child nodes


    def name(self):
        """The tag name of the node."""
        return self._name

    def nChildren(self):
        """Number of child elements."""
        return len(self._children)

    def addChild(self, name, value=""):
        """Add a child with tag 'name', and set its value if the value
        parameter is supplied."""

        # create a new node for the child
        c = XMLnode(name = name, value = value)

        # add it to the list of children, and to the dictionary
        # of children
        self._children.append(c)
        self._childmap[name] = c
        return c

    def addComment(self, comment):
        """Add a comment."""
        self.addChild(name = '_comment_', value = comment)

    def value(self):
        """A string containing the element value."""
        return self._value

    def child(self, name=""):
        """The child node with specified name."""
        return self._childmap[name]

    def children(self):
        """ An iterator over the child nodes """
        for c in self._children:
            yield c

    def __getitem__(self, key):
        """Get an attribute using the syntax node[key]"""
        return self._attribs[key]

    def __setitem__(self, key, value):
        """Set a new attribute using the syntax node[key] = value."""
        self._attribs[key] = value

    def __call__(self):
        """Allows getting the value using the syntax 'node()'"""
        return self._value

    def write(self, filename):
        """Write out the XML tree to a file."""
        s = ['<?xml version="1.0"?>\n']
        self._write(s, 0)
        s.append('\n')
        with open(filename, 'w') as f:
            f.write(''.join(s))

    def write_comment(self, s, level):
        s.append('\n'+indent[level]+'<!--')
        value = self._value
        if value:
            if value[0] != ' ':
                value = ' '+value
            if value[-1] != ' ':
                value = value+' '
        s.append(value+'-->')

    def write_attribs(self, s):
        for a in self._attribs:
            s.append(' '+a+'="'+self._attribs[a]+'"')

    def write_value(self, s, level):
        indnt = indent[level]
        vv = self._value.lstrip()
        ieol = vv.find('\n')
        if ieol >= 0:
            while True:
                ieol = vv.find('\n')
                if ieol >= 0:
                    s.extend(('\n  ', indnt, vv[:ieol]))
                    vv = vv[ieol+1:].lstrip()
                else:
                    s.extend(('\n  ',indnt,vv))
                    break
        else:
            s.append(self._value)

    def _write(self, s, level = 0):
        """Internal method used to write the XML representation of each node."""
        if not self.name:
            return

        # handle comments
        if self._name == '_comment_':
            self.write_comment(s, level)
            return

        indnt = indent[level]

        # write the opening tag and attributes
        s.extend((indnt, '<', self._name))
        self.write_attribs(s)

        if not self._value and not self._children:
            s.append('/>')
        else:
            s.append('>')
            if self._value:
                self.write_value(s, level)

            for c in self._children:
                s.append('\n')
                c._write(s, level + 2)
            if self._children:
                s.extend(('\n', indnt))
            s.extend(('</', self._name, '>'))

#--------------------------------------------------

# constants that can be used in .cti files
OneAtm = 1.01325e5
OneBar = 1.0e5
# Conversion from eV to J/kmol (electronCharge * Navrog)
eV = 9.64853364595687e7
# Electron Mass in kg
ElectronMass = 9.10938291e-31

import math, copy

# default units
_ulen = 'm'
_umol = 'kmol'
_umass = 'kg'
_utime = 's'
_ue = 'J/kmol'
_uenergy = 'J'
_upres = 'Pa'

# used to convert reaction pre-exponentials
_length = {'cm':0.01, 'm':1.0, 'mm':0.001}
_moles = {'kmol':1.0, 'mol':0.001, 'molec':1.0/6.02214129e26}
_time = {'s':1.0, 'min':60.0, 'hr':3600.0}

# default std state pressure
_pref = 1.0e5    # 1 bar

_name = 'noname'

# these lists store top-level entries
_elements = []
_species = []
_speciesnames = []
_phases = []
_reactions = []
_atw = {}
_enames = {}

_valsp = ''
_valrxn = ''
_valexport = ''
_valfmt = ''

def export_species(filename, fmt = 'CSV'):
    global _valexport
    global _valfmt
    _valexport = filename
    _valfmt = fmt

def validate(species = 'yes', reactions = 'yes'):
    global _valsp
    global _valrxn
    _valsp = species
    _valrxn = reactions

def isnum(a):
    """True if a is an integer or floating-point number."""
    if isinstance(a, (int, float)):
        return 1
    else:
        return 0

def is_local_species(name):
    """true if the species named 'name' is defined in this file"""
    if name in _speciesnames:
        return 1
    return 0

def dataset(nm):
    "Set the dataset name. Invoke this to change the name of the xml file."
    global _name
    _name = nm

def standard_pressure(p0):
    """Set the default standard-state pressure."""
    global _pref
    _pref = p0

def units(length = '', quantity = '', mass = '', time = '',
          act_energy = '', energy = '', pressure = ''):
    """
    Set the default units.

    :param length:
        The default units for length. Default: ``'m'``
    :param mass:
        The default units for mass. Default: ``'kg'``
    :param quantity:
        The default units to specify number of molecules. Default: ``'kmol'``
    :param time:
        The default units for time. Default: ``'s'``
    :param energy:
        The default units for energies. Default: ``'J'``
    :param act_energy:
        The default units for activation energies. Default: ``'K'``
    :param pressure:
        The default units for pressure. Default: ``'Pa'``
    """
    global _ulen, _umol, _ue, _utime, _umass, _uenergy, _upres
    if length: _ulen = length
    if quantity: _umol = quantity
    if act_energy: _ue = act_energy
    if time: _utime = time
    if mass: _umass = mass
    if energy: _uenergy = energy
    if pressure: _upres = pressure

def ufmt(base, n):
    """return a string representing a unit to a power n."""
    if n == 0: return ''
    if n == 1: return '-'+base
    if n == -1: return '/'+base
    if n > 0: return '-'+base+str(n)
    if n < 0: return '/'+base+str(-n)

def write(outName=None):
    """write the CTML file."""
    x = XMLnode("ctml")
    v = x.addChild("validate")
    v["species"] = _valsp
    v["reactions"] = _valrxn

    if _elements:
        ed = x.addChild("elementData")
        for e in _elements:
            e.build(ed)

    for ph in _phases:
        ph.build(x)
    s = species_set(name = _name, species = _species)
    s.build(x)

    r = x.addChild('reactionData')
    r['id'] = 'reaction_data'
    for rx in _reactions:
        rx.build(r)

    if outName is not None:
        x.write(outName)
    elif _name != 'noname':
        x.write(_name+'.xml')
    else:
        print(x)

    if _valexport:
        f = open(_valexport,'w')
        for s in _species:
            s.export(f, _valfmt)
        f.close()

def addFloat(x, nm, val, fmt='', defunits=''):
    """
    Add a child element to XML element x representing a
    floating-point number.
    """
    u = ''
    s = ''
    if isnum(val):
        fval = float(val)
        if fmt:
            s = fmt % fval
        else:
            s = repr(fval)
        xc = x.addChild(nm, s)
        if defunits:
            xc['units'] = defunits
    else:
        v = val[0]
        u = val[1]
        if fmt:
            s = fmt % v
        else:
            s = repr(v)
        xc = x.addChild(nm, s)
        xc['units'] = u


def getAtomicComp(atoms):
    if isinstance(atoms, dict): return atoms
    a = atoms.replace(',',' ')
    toks = a.split()
    d = {}
    for t in toks:
        b = t.split(':')
        d[b[0]] = int(b[1])
    return d

def getReactionSpecies(s):
    """Take a reaction string and return a
    dictionary mapping species names to stoichiometric
    coefficients. If any species appears more than once,
    the returned stoichiometric coefficient is the sum.
    >>> s = 'CH3 + 3 H + 5.2 O2 + 0.7 H'
    >>> getReactionSpecies(s)
    >>> {'CH3':1, 'H':3.7, 'O2':5.2}
    """

    # get rid of the '+' signs separating species. Only plus signs
    # surrounded by spaces are replaced, so that plus signs may be
    # used in species names (e.g. 'Ar3+')
    toks = s.replace(' + ',' ').split()
    d = {}
    n = 1.0
    for t in toks:

        # try to convert the token to a number.
        try:
            n = float(t)
            if n < 0.0:
                raise CTI_Error("negative stoichiometric coefficient:"
                                +s)

        #if t > '0' and t < '9':
        #    n = int(t)
        #else:

        # token isn't a number, so it must be a species name
        except:
            # already seen this token so increment its value by the last
            # value of n
            if t in d:
                d[t] += n
            else:
                # first time this token has been seen, so set its value to n
                d[t] = n

            # reset n to 1.0 for species that do not specify a stoichiometric
            # coefficient
            n = 1

    return d


class element(object):
    """ An atomic element or isotope. """
    def __init__(self, symbol = '',
                 atomic_mass = 0.01,
                 atomic_number = 0):
        """
        :param symbol:
            The symbol for the element or isotope.
        :param atomic_mass:
            The atomic mass in amu.
        """
        self._sym = symbol
        self._atw = atomic_mass
        self._num = atomic_number
        global _elements
        _elements.append(self)

    def build(self, db):
        e = db.addChild("element")
        e["name"] = self._sym
        e["atomicWt"] = repr(self._atw)
        e["atomicNumber"] = repr(self._num)


class species_set(object):
    def __init__(self, name = '', species = []):
        self._s = species
        self._name = name
        #self.type = SPECIES_SET

    def build(self, p):
        p.addComment('     species definitions     ')
        sd = p.addChild("speciesData")
        sd["id"] = "species_data"
        for s in self._s:
            #if s.type == SPECIES:
            s.build(sd)
            #else:
            #    raise 'wrong object type in species_set: '+s.__class__


class species(object):
    """A constituent of a phase or interface."""

    def __init__(self,
                 name = 'missing name!',
                 atoms = '',
                 note = '',
                 thermo = None,
                 transport = None,
                 charge = -999,
                 size = 1.0):
        """
        :param name:
            The species name (or formula). The name may be arbitrarily long,
            although usually a relatively short, abbreviated name is most
            convenient. Required parameter.
        :param atoms:
            The atomic composition, specified by a string containing
            space-delimited <element>:<atoms> pairs. The number of atoms may be
            either an integer or a floating-point number.
        :param note:
            A user-defined comment. Not evaluated by Cantera itself.
        :param thermo:
            The parameterization to use to compute the reference-state
            thermodynamic properties. This must be one of the entry types
            described in :ref:`sec-thermo-models`. To specify multiple
            parameterizations, each for a different temperature range,
            group them in parentheses.
        :param transport:
            An entry specifying parameters to compute this species'
            contribution to the transport properties. This must be one of the
            entry types described in :ref:`sec-species-transport-models`, and
            must be consistent with the transport model of the phase into which
            the species is imported. To specify parameters for multiple
            transport models, group the entries in parentheses.
        :param size:
            The species "size". Currently used only for surface species,
            where it represents the number of sites occupied.
        :param charge:
            The charge, in multiples of :math:`|e|`. If not specified, the
            charge will be calculated from the number of "atoms" of element
            ``E``, which represents an electron.
        """
        self._name = name
        self._atoms = getAtomicComp(atoms)

        self._comment = note

        if thermo:
            self._thermo = thermo
        else:
            self._thermo = const_cp()

        self._transport = transport
        chrg = 0
        self._charge = charge
        if 'E' in self._atoms:
            chrg = -self._atoms['E']
            if self._charge != -999:
                if self._charge != chrg:
                    raise CTI_Error('specified charge inconsistent with number of electrons')
            else:
                self._charge = chrg
        self._size = size

        global _species
        global _enames
        _species.append(self)
        global _speciesnames
        if name in _speciesnames:
            raise CTI_Error('species '+name+' multiply defined.')
        _speciesnames.append(name)
        for e in self._atoms.keys():
            _enames[e] = 1

    def export(self, f, fmt = 'CSV'):
        global _enames
        if fmt == 'CSV':
            s = self._name+','
            for e in _enames:
                if e in self._atoms:
                    s += repr(self._atoms[e])+','
                else:
                    s += '0,'
            f.write(s)
            if isinstance(self._thermo, thermo):
                self._thermo.export(f, fmt)
            else:
                nt = len(self._thermo)
                for n in range(nt):
                    self._thermo[n].export(f, fmt)
            f.write('\n')


    def build(self, p):
        hdr = '    species '+self._name+'    '
        p.addComment(hdr)
        s = p.addChild("species")
        s["name"] = self._name
        a = ''
        for e in self._atoms.keys():
            a += e+':'+str(self._atoms[e])+' '
        s.addChild("atomArray",a)
        if self._comment:
            s.addChild("note",self._comment)
        if self._charge != -999:
            s.addChild("charge",self._charge)
        if self._size != 1.0:
            s.addChild("size",self._size)
        if self._thermo:
            t = s.addChild("thermo")
            if isinstance(self._thermo, thermo):
                self._thermo.build(t)
            else:
                nt = len(self._thermo)
                for n in range(nt):
                    self._thermo[n].build(t)
        if self._transport:
            t = s.addChild("transport")
            if isinstance(self._transport, transport):
                self._transport.build(t)
            else:
                nt = len(self._transport)
                for n in range(nt):
                    self._transport[n].build(t)

class thermo(object):
    """Base class for species standard-state thermodynamic properties."""
    def _build(self, p):
        return p.addChild("thermo")
    def export(self, f, fmt = 'CSV'):
        pass

class Mu0_table(thermo):
    """Properties are computed by specifying a table of standard
    chemical potentials vs. T."""

    def __init__(self, Trange = (0.0, 0.0),
                 h298 = 0.0,
                 mu0 = None,
                 p0 = -1.0):
        self._t = Trange
        self._h298 = h298
        self._mu0 = mu0
        self._pref = p0

    def build(self, t):
        n = t.addChild("Mu0")
        n['Tmin'] = repr(self._t[0])
        n['Tmax'] = repr(self._t[1])
        if self._pref <= 0.0:
            n['P0'] = repr(_pref)
        else:
            n['P0'] = repr(self._pref)
        energy_units = _uenergy+'/'+_umol
        addFloat(n,"H298", self._h298, defunits = energy_units)
        n.addChild("numPoints", len(self._mu0))


        mustr = ''
        tstr = ''
        col = 0
        for v in self._mu0:
            mu0 = v[1]
            t = v[0]
            tstr += '%17.9E, ' % t
            mustr += '%17.9E, ' % mu0
            col += 1
            if col == 3:
                tstr = tstr[:-2]+'\n'
                mustr = mustr[:-2]+'\n'
                col = 0

        u = n.addChild("floatArray", mustr)
        u["size"] = "numPoints"
        u["name"] = "Mu0Values"

        u = n.addChild("floatArray", tstr)
        u["size"] = "numPoints"
        u["name"] = "Mu0Temperatures"


class NASA(thermo):
    """The 7-coefficient NASA polynomial parameterization."""
    def __init__(self, Trange = (0.0, 0.0), coeffs = [], p0 = -1.0):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required.
        :param coeffs:
            Array of seven coefficients :math:`(a_0, \ldots , a_6)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value is used, which is set by the ``standard_pressure``
            directive.
        """
        self._t = Trange
        self._pref = p0
        if len(coeffs) != 7:
            raise CTI_Error('NASA coefficient list must have length = 7')
        self._coeffs = coeffs


    def export(self, f, fmt='CSV'):
        if fmt == 'CSV':
            s = 'NASA,'+str(self._t[0])+','+str(self._t[1])+','
            for i in  range(7):
                s += '%17.9E, ' % self._coeffs[i]
            f.write(s)

    def build(self, t):
        n = t.addChild("NASA")
        n['Tmin'] = repr(self._t[0])
        #n['Tmid'] = repr(self._t[1])
        n['Tmax'] = repr(self._t[1])
        if self._pref <= 0.0:
            n['P0'] = repr(_pref)
        else:
            n['P0'] = repr(self._pref)
        s = ''
        for i in  range(4):
            s += '%17.9E, ' % self._coeffs[i]
        s += '\n'
        s += '%17.9E, %17.9E, %17.9E' % (self._coeffs[4],
                                           self._coeffs[5], self._coeffs[6])
        #if i > 0 and 3*((i+1)/3) == i: s += '\n'
        #s = s[:-2]
        u = n.addChild("floatArray", s)
        u["size"] = "7"
        u["name"] = "coeffs"


class NASA9(thermo):
    """NASA9 polynomial parameterization for a single temperature region."""

    def __init__(self, Trange = (0.0, 0.0),
                 coeffs = [], p0 = -1.0):
        self._t = Trange         # Range of the polynomial representation
        self._pref = p0          # Reference pressure
        if len(coeffs) != 9:
            raise CTI_Error('NASA9 coefficient list must have length = 9')
        self._coeffs = coeffs


    def export(self, f, fmt='CSV'):
        if fmt == 'CSV':
            s = 'NASA9,'+str(self._t[0])+','+str(self._t[1])+','
            for i in  range(9):
                s += '%17.9E, ' % self._coeffs[i]
            f.write(s)

    def build(self, t):
        n = t.addChild("NASA9")
        n['Tmin'] = repr(self._t[0])
        n['Tmax'] = repr(self._t[1])
        if self._pref <= 0.0:
            n['P0'] = repr(_pref)
        else:
            n['P0'] = repr(self._pref)
        s = ''
        for i in  range(4):
            s += '%17.9E, ' % self._coeffs[i]
        s += '\n'
        s += '%17.9E, %17.9E, %17.9E, %17.9E,' % (self._coeffs[4], self._coeffs[5],
                                                    self._coeffs[6], self._coeffs[7])
        s += '\n'
        s += '%17.9E' % (self._coeffs[8])
        u = n.addChild("floatArray", s)
        u["size"] = "9"
        u["name"] = "coeffs"


class activityCoefficients(object):
    pass


class pureFluidParameters(activityCoefficients):
    """
    """

    def __init__(self, species = None, a_coeff = [], b_coeff = 0):
        """
        """
        self._species = species
        self._acoeff = a_coeff
        self._bcoeff = b_coeff

    def build(self,a):
        f= a.addChild("pureFluidParameters")
        f['species'] = self._species
        s = '%10.4E, %10.4E \n' % (self._acoeff[0], self._acoeff[1])
        ac = f.addChild("a_coeff",s)
        ac["units"] = _upres+'-'+_ulen+'6/'+_umol+'2'
        ac["model"] = "linear_a"
        s = '%0.2f \n' % self._bcoeff
        bc = f.addChild("b_coeff",s)
        bc["units"] = _ulen+'3/'+_umol


class crossFluidParameters(activityCoefficients):
    def __init__(self, species = None, a_coeff = [], b_coeff = []):
        self._species1, self._species2 = species.split(' ')
        self._acoeff = a_coeff
        self._bcoeff = b_coeff

    def build(self,a):
        f= a.addChild("crossFluidParameters")
        f["species2"] = self._species2
        f["species1"] = self._species1
        s = '%10.4E, %10.4E \n' % (self._acoeff[0], self._acoeff[1])
        ac = f.addChild("a_coeff",s)
        ac["units"] = _upres+'-'+_ulen+'6/'+_umol+'2'
        ac["model"] = "linear_a"
        if self._bcoeff:
            s = '%0.2f \n' % self._bcoeff
            bc = f.addChild("b_coeff",s)
            bc["units"] = _ulen+'3/'+_umol


class Shomate(thermo):
    """Shomate polynomial parameterization."""

    def __init__(self, Trange = (0.0, 0.0), coeffs = [], p0 = -1.0):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required input.
        :param coeffs:
            Sequence of seven coefficients :math:`(A, \ldots ,G)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value set by the ``standard_pressure`` directive is used.
        """
        self._t = Trange
        self._pref = p0
        if len(coeffs) != 7:
            raise CTI_Error('Shomate coefficient list must have length = 7')
        self._coeffs = coeffs


    def build(self, t):
        n = t.addChild("Shomate")
        n['Tmin'] = repr(self._t[0])
        n['Tmax'] = repr(self._t[1])
        if self._pref <= 0.0:
            n['P0'] = repr(_pref)
        else:
            n['P0'] = repr(self._pref)
        s = ''
        for i in  range(4):
            s += '%17.9E, ' % self._coeffs[i]
        s += '\n'
        s += '%17.9E, %17.9E, %17.9E' % (self._coeffs[4],
                                         self._coeffs[5], self._coeffs[6])
        u = n.addChild("floatArray", s)
        u["size"] = "7"
        u["name"] = "coeffs"


class Adsorbate(thermo):
    """Adsorbed species characterized by a binding energy and a set of
    vibrational frequencies."""

    def __init__(self, Trange = (0.0, 0.0),
                 binding_energy = 0.0,
                 frequencies = [], p0 = -1.0):
        self._t = Trange
        self._pref = p0
        self._freqs = frequencies
        self._be = binding_energy


    def build(self, t):
        n = t.addChild("adsorbate")
        n['Tmin'] = repr(self._t[0])
        n['Tmax'] = repr(self._t[1])
        if self._pref <= 0.0:
            n['P0'] = repr(_pref)
        else:
            n['P0'] = repr(self._pref)

        energy_units = _uenergy+'/'+_umol
        addFloat(n,'binding_energy',self._be, defunits = energy_units)
        s = ""
        nfreq = len(self._freqs)
        for i in  range(nfreq):
            s += '%17.9E, ' % self._freqs[i]
        s += '\n'
        u = n.addChild("floatArray", s)
        u["size"] = repr(nfreq)
        u["name"] = "freqs"



class const_cp(thermo):
    """Constant specific heat."""

    def __init__(self,
                 t0 = 298.15, cp0 = 0.0, h0 = 0.0, s0 = 0.0,
                 tmax = 5000.0, tmin = 100.0):
        """
        :param t0:
            Temperature parameter T0. Default: 298.15 K.
        :param cp0:
            Reference-state molar heat capacity (constant). Default: 0.0.
        :param h0:
            Reference-state molar enthalpy at temperature T0. Default: 0.0.
        :param s0:
            Reference-state molar entropy at temperature T0. Default: 0.0.
        """
        self._t = [tmin, tmax]
        self._c = [t0, h0, s0, cp0]

    def build(self, t):
        #t = self._build(p)
        c = t.addChild('const_cp')
        if self._t[0] >= 0.0: c['Tmin'] = repr(self._t[0])
        if self._t[1] >= 0.0: c['Tmax'] = repr(self._t[1])
        energy_units = _uenergy+'/'+_umol
        addFloat(c,'t0',self._c[0], defunits = 'K')
        addFloat(c,'h0',self._c[1], defunits = energy_units)
        addFloat(c,'s0',self._c[2], defunits = energy_units+'/K')
        addFloat(c,'cp0',self._c[3], defunits = energy_units+'/K')

class transport(object):
    pass

class gas_transport(transport):
    """
    Species-specific Transport coefficients for ideal gas transport models.
    """
    def __init__(self, geom = 'nonlin',
                 diam = 0.0, well_depth = 0.0, dipole = 0.0,
                 polar = 0.0, rot_relax = 0.0):
        """
        :param geom:
            A string specifying the molecular geometry. One of ``atom``,
            ``linear``, or ``nonlin``. Required.
        :param diam:
            The Lennard-Jones collision diameter in Angstroms. Required.
        :param well_depth:
            The Lennard-Jones well depth in Kelvin. Required.
        :param dipole:
            The permanent dipole moment in Debye. Default: 0.0
        :param polar:
            The polarizability in A^3. Default: 0.0
        :param rot_relax:
            The rotational relaxation collision number at 298 K. Dimensionless.
            Default: 0.0
        """
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

class rate_expression(object):
    pass

class Arrhenius(rate_expression):
    def __init__(self,
                 A = 0.0,
                 n = 0.0,
                 E = 0.0,
                 coverage = [],
                 rate_type = ''):
        """
        :param A:
            The pre-exponential coefficient. Required input. If entered without
            units, the units will be computed considering all factors that
            affect the units. The resulting units string is written to the CTML
            file individually for each reaction pre-exponential coefficient.
        :param n:
            The temperature exponent. Dimensionless. Default: 0.0.
        :param E:
            Activation energy. Default: 0.0.
        """
        self._c = [A, n, E]
        self._type = rate_type

        if coverage:
            if isinstance(coverage[0], str):
                self._cov = [coverage]
            else:
                self._cov = coverage
        else:
            self._cov = None


    def build(self, p, units_factor = 1.0,
              gas_species = [], name = '', rxn_phase = None):

        a = p.addChild('Arrhenius')
        if name: a['name'] = name

        # check for sticking probability
        if self._type:
            a['type'] = self._type
            if self._type == 'stick':
                ngas = len(gas_species)
                if ngas != 1:
                    raise CTI_Error("""
Sticking probabilities can only be used for reactions with one gas-phase
reactant, but this reaction has """+str(ngas)+': '+str(gas_species))
                else:
                    a['species'] = gas_species[0]
                    units_factor = 1.0


        # if a pure number is entered for A, multiply by the conversion
        # factor to SI and write it to CTML as a pure number. Otherwise,
        # pass it as-is through to CTML with the unit string.
        if isnum(self._c[0]):
            addFloat(a,'A',self._c[0]*units_factor, fmt = '%14.6E')
        elif len(self._c[0]) == 2 and self._c[0][1] == '/site':
            addFloat(a,'A',self._c[0][0]/rxn_phase._sitedens,
                     fmt = '%14.6E')
        else:
            addFloat(a,'A',self._c[0], fmt = '%14.6E')

        # The b coefficient should be dimensionless, so there is no
        # need to use 'addFloat'
        a.addChild('b', repr(self._c[1]))

        # If a pure number is entered for the activation energy,
        # add the default units, otherwise use the supplied units.
        addFloat(a,'E', self._c[2], fmt = '%f', defunits = _ue)

        # for surface reactions, a coverage dependence may be specified.
        if self._cov:
            for cov in self._cov:
                c = a.addChild('coverage')
                c['species'] = cov[0]
                addFloat(c, 'a', cov[1], fmt = '%f')
                c.addChild('m', repr(cov[2]))
                addFloat(c, 'e', cov[3], fmt = '%f', defunits = _ue)

def stick(A = 0.0, n = 0.0, E = 0.0, coverage = []):
    return Arrhenius(A = A, n = n, E = E, coverage = coverage, rate_type = 'stick')


def getPairs(s):
    toks = s.split()
    m = {}
    for t in toks:
        key, val = t.split(':')
        m[key] = float(val)
    return m

class reaction(object):
    """
    A homogeneous chemical reaction with pressure-independent rate coefficient
    and mass-action kinetics.
    """
    def __init__(self,
                 equation = '',
                 kf = None,
                 id = '',
                 order = '',
                 options = []):
        """
        :param equation:
            A string specifying the chemical equation.
        :param rate_coeff:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, n,E] in
            the modified Arrhenius function :math:`A T^n exp(-E/\hat{R}T)`.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options:
            Processing options, as described in :ref:`sec-phase-options`.
        """
        self._id = id
        self._e = equation
        self._order = order

        if isinstance(options, str):
            self._options = [options]
        else:
            self._options = options
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

        self._rxnorder = copy.copy(self._r)
        if self._order:
            order = getPairs(self._order)
            for o in order.keys():
                if o in self._rxnorder:
                    self._rxnorder[o] = order[o]
                else:
                    raise CTI_Error("order specified for non-reactant: "+o)

        self._kf = kf
        self._igspecies = []
        self._dims = [0]*4
        self._rxnphase = None
        self._type = ''
        _reactions.append(self)

    def unit_factor(self):
        """
        Conversion factor from given rate constant units to the MKS (+kmol)
        used internally by Cantera, taking into account the reaction order.
        """
        return (math.pow(_length[_ulen], -self.ldim) *
                math.pow(_moles[_umol], -self.mdim) / _time[_utime])

    def build(self, p):
        if self._id:
            id = self._id
        else:
            id = '%04i' % self._num

        self.mdim = 0
        self.ldim = 0

        rxnph = []
        for s in self._r:
            ns = self._rxnorder[s]
            nm = -999
            nl = -999

            mindim = 4
            for ph in _phases:
                if ph.has_species(s):
                    nm, nl = ph.conc_dim()
                    if ph.is_ideal_gas():
                        self._igspecies.append(s)
                    if not ph in rxnph:
                        rxnph.append(ph)
                        self._dims[ph._dim] += 1
                        if ph._dim < mindim:
                            self._rxnphase = ph
                            mindim = ph._dim
                    break
            if nm == -999:
                raise CTI_Error("species "+s+" not found")

            self.mdim += nm*ns
            self.ldim += nl*ns

        p.addComment("   reaction "+id+"    ")
        r = p.addChild('reaction')
        r['id'] = id
        if self.rev:
            r['reversible'] = 'yes'
        else:
            r['reversible'] = 'no'

        for s in self._options:
            if s == 'duplicate':
                r['duplicate'] = 'yes'
            elif s == 'negative_A':
                r['negative_A'] = 'yes'

        ee = self._e.replace('<','[').replace('>',']')
        r.addChild('equation',ee)

        if self._order:
            for osp in self._rxnorder:
                o = r.addChild('order',self._rxnorder[osp])
                o['species'] = osp


        # adjust the moles and length powers based on the dimensions of
        # the rate of progress (moles/length^2 or moles/length^3)
        if self._type == 'surface':
            self.mdim += -1
            self.ldim += 2
            p = self._dims[:3]
            if p[0] != 0 or p[1] != 0 or p[2] > 1:
                raise CTI_Error(self._e +'\nA surface reaction may contain at most '+
                                   'one surface phase.')
        elif self._type == 'edge':
            self.mdim += -1
            self.ldim += 1
            p = self._dims[:2]
            if p[0] != 0 or p[1] > 1:
                raise CTI_Error(self._e+'\nAn edge reaction may contain at most '+
                                   'one edge phase.')
        else:
            self.mdim += -1
            self.ldim += 3

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
        elif self._type == 'edge':
            self._kf = [self._kf]
        elif self._type == 'threeBody':
            self._kf = [self._kf]
            self.mdim += 1
            self.ldim -= 3
        elif self._type == 'chebyshev':
            self._kf = []

        if self._type == 'edge':
            if self._beta > 0:
                electro = kfnode.addChild('electrochem')
                electro['beta'] = repr(self._beta)

        for kf in self._kf:
            if isinstance(kf, rate_expression):
                k = kf
            else:
                k = Arrhenius(A = kf[0], n = kf[1], E = kf[2])
            k.build(kfnode, self.unit_factor(), gas_species = self._igspecies,
                    name = nm, rxn_phase = self._rxnphase)

            if self._type == 'falloff':
                # set values for low-pressure rate coeff if falloff rxn
                self.mdim += 1
                self.ldim -= 3
                nm = 'k0'
            elif self._type == 'chemAct':
                # set values for high-pressure rate coeff if this is a
                # chemically activated reaction
                self.mdim -= 1
                self.ldim += 3
                nm = 'kHigh'

        rstr = ' '.join('%s:%s' % item for item in self._r.items())
        pstr = ' '.join('%s:%s' % item for item in self._p.items())
        r.addChild('reactants',rstr)
        r.addChild('products', pstr)
        return r

#-------------------


class three_body_reaction(reaction):
    """
    A three-body reaction.
    """
    def __init__(self,
                 equation = '',
                 kf = None,
                 efficiencies = '',
                 id = '',
                 options = []
                 ):
        """
        :param equation:
            A string specifying the chemical equation. The reaction can be
            written in either the association or dissociation directions, and
            may be reversible or irreversible.
        :param rate_coeff:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A,n,E] in
            the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies.
            The efficiencies for unspecified species are set to 1.0.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options:
            Processing options, as described in :ref:`sec-phase-options`.
        """
        reaction.__init__(self, equation, kf, id, '', options)
        self._type = 'threeBody'
        self._effm = 1.0
        self._eff = efficiencies

        # clean up reactant and product lists
        for r in list(self._r.keys()):
            if r == 'M' or r == 'm':
                del self._r[r]
        for p in list(self._p.keys()):
            if p == 'M' or p == 'm':
                del self._p[p]

    def build(self, p):
        r = reaction.build(self, p)
        if r == 0: return
        kfnode = r.child('rateCoeff')

        if self._eff:
            eff = kfnode.addChild('efficiencies',self._eff)
            eff['default'] = repr(self._effm)


class pdep_reaction(reaction):
    """ Base class for falloff_reaction and chemically_activated_reaction """

    def clean_up_reactants_products(self):
        del self._r['(+']
        del self._p['(+']
        if 'M)' in self._r:
            del self._r['M)']
            del self._p['M)']
        elif 'm)' in self._r:
            del self._r['m)']
            del self._p['m)']
        else:
            for r in list(self._r.keys()):
                if r[-1] == ')' and r.find('(') < 0:
                    species = r[:-1]
                    if self._eff:
                        raise CTI_Error('(+ '+species+') and '+self._eff+' cannot both be specified')
                    self._eff = species+':1.0'
                    self._effm = 0.0

                    del self._r[r]
                    del self._p[r]

    def build(self, p):
        r = reaction.build(self, p)
        if r == 0: return
        kfnode = r.child('rateCoeff')

        if self._eff and self._effm >= 0.0:
            eff = kfnode.addChild('efficiencies',self._eff)
            eff['default'] = repr(self._effm)

        if self._falloff:
            self._falloff.build(kfnode)


class falloff_reaction(pdep_reaction):
    """ A gas-phase falloff reaction. """
    def __init__(self, equation, kf0, kf,
                 efficiencies='', falloff=None, id='', options=[]):
        """
        :param equation:
            A string specifying the chemical equation.
        :param rate_coeff_inf:
            The rate coefficient for the forward direction in the high-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, n,E] in the modified Arrhenius function.
        :param rate_coeff_0:
            The rate coefficient for the forward direction in the low-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, n,E] in the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies. The
            efficiency for unspecified species is set to 1.0.
        :param falloff:
            An embedded entry specifying a falloff function. If omitted, a
            unity falloff function (Lindemann form) will be used.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options:
            Processing options, as described in :ref:`sec-phase-options`.
        """
        kf2 = (kf, kf0)
        reaction.__init__(self, equation, kf2, id, '', options)
        self._type = 'falloff'
        # use a Lindemann falloff function by default
        self._falloff = falloff
        if self._falloff == None:
            self._falloff = Lindemann()

        self._effm = 1.0
        self._eff = efficiencies

        self.clean_up_reactants_products()


class chemically_activated_reaction(pdep_reaction):
    """ A gas-phase, chemically activated reaction. """

    def __init__(self, equation, kLow, kHigh,
                 efficiencies='', falloff=None, id='', options=[]):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kLow:
            The rate coefficient for the forward direction in the low-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, n,E] in the modified Arrhenius function.
        :param kHigh:
            The rate coefficient for the forward direction in the high-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, n,E] in the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies. The
            efficiency for unspecified species is set to 1.0.
        :param falloff:
            An embedded entry specifying a falloff function. If omitted, a
            unity falloff function (Lindemann form) will be used.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options:
            Processing options, as described in :ref:`sec-phase-options`.
        """
        reaction.__init__(self, equation, (kLow, kHigh), id, '', options)
        self._type = 'chemAct'
        # use a Lindemann falloff function by default
        self._falloff = falloff
        if self._falloff == None:
            self._falloff = Lindemann()

        self._effm = 1.0
        self._eff = efficiencies

        self.clean_up_reactants_products()

class pdep_arrhenius(reaction):
    """
    Pressure-dependent rate calculated by interpolating between Arrhenius
    expressions at different pressures.

    :param equation:
        A string specifying the chemical equation.
    :param args:
        Each additional argument is a sequence of four elements specifying the
        pressure and the Arrhenius parameters at that pressure.
    """
    def __init__(self, equation='', *args, **kwargs):
        self.pressures = []
        self.arrhenius = []
        for p, A, n, Ea in args:
            self.pressures.append(p)
            self.arrhenius.append((A, n, Ea))

        reaction.__init__(self, equation, self.arrhenius, **kwargs)
        self._type = 'plog'

    def build(self, p):
        r = reaction.build(self, p)
        kfnode = r.child('rateCoeff')
        for i,c in enumerate(kfnode.children()):
            assert c.name() == 'Arrhenius'
            addFloat(c, 'P', self.pressures[i])


class chebyshev_reaction(reaction):
    """
    Pressure-dependent rate calculated in terms of a bivariate Chebyshev
    polynomial.

    :param equation:
        A string specifying the chemical equation.
    :param Tmin:
        The minimum temperature at which the rate expression is defined
    :param Tmax:
        the maximum temperature at which the rate expression is defined
    :param Pmin:
        The minimum pressure at which the rate expression is defined
    :param Pmax:
        The maximum pressure at which the rate expression is defined
    :param coeffs:
        A 2D array of the coefficients defining the rate expression. For a
        polynomial with M points in temperature and N points in pressure, this
        should be a list of M lists each with N elements.
    """
    def __init__(self, equation='', Tmin=300.0, Tmax=2500.0,
                 Pmin=(0.001, 'atm'), Pmax=(100.0, 'atm'),
                 coeffs=[[]], **kwargs):
        reaction.__init__(self, equation, **kwargs)
        self._type = 'chebyshev'
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.coeffs = coeffs

        # clean up reactant and product lists
        del self._r['(+']
        del self._p['(+']
        if 'M)' in self._r:
            del self._r['M)']
            del self._p['M)']
        if 'm)' in self._r:
            del self._r['m)']
            del self._p['m)']

    def build(self, p):
        r = reaction.build(self, p)
        kfnode = r.child('rateCoeff')
        addFloat(kfnode, 'Tmin', self.Tmin)
        addFloat(kfnode, 'Tmax', self.Tmax)
        addFloat(kfnode, 'Pmin', self.Pmin)
        addFloat(kfnode, 'Pmax', self.Pmax)

        self.coeffs[0][0] += math.log10(self.unit_factor());

        lines = []
        for line in self.coeffs:
            lines.append(', '.join('{0:12.5e}'.format(val)
                                   for val in line))

        coeffNode = kfnode.addChild('floatArray', ',\n'.join(lines))
        coeffNode['name'] = 'coeffs'
        coeffNode['degreeT'] = str(len(self.coeffs))
        coeffNode['degreeP'] = str(len(self.coeffs[0]))


class surface_reaction(reaction):
    """
    A heterogeneous chemical reaction with pressure-independent rate
    coefficient and mass-action kinetics.
    """
    def __init__(self, equation='', kf=None, id='', order='', options=[]):
        """
        :param equation:
            A string specifying the chemical equation.
        :param rate_coeff:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, n,E] in
            the modified Arrhenius function.
        :param sticking_prob:
            The reactive sticking probability for the forward direction. This
            can only be specified if there is only one bulk-phase reactant and
            it belongs to an ideal gas phase. If a sequence of three numbers is
            given, these will be interpreted as [A, n,E] in the modified
            Arrhenius function.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options:
            Processing options, as described in :ref:`sec-phase-options`.
        """
        reaction.__init__(self, equation, kf, id, order, options)
        self._type = 'surface'


class edge_reaction(reaction):

    def __init__(self,
                 equation = '',
                 kf = None,
                 id = '',
                 order = '',
                 beta = 0.0,
                 options = []):
        reaction.__init__(self, equation, kf, id, order, options)
        self._type = 'edge'
        self._beta = beta


#--------------


class state(object):
    """
    An embedded entry that specifies the thermodynamic state of a phase
    or interface.
    """
    def __init__(self,
                 temperature = None,
                 pressure = None,
                 mole_fractions = None,
                 mass_fractions = None,
                 density = None,
                 coverages = None,
                 solute_molalities = None):
        """
        :param temperature:
            The temperature.
        :param pressure:
            The pressure.
        :param density:
            The density. Cannot be specified if the phase is incompressible.
        :param mole_fractions:
            A string specifying the species mole fractions. Unspecified species
            are set to zero.
        :param mass_fractions:
            A string specifying the species mass fractions. Unspecified species
            are set to zero.
        :param coverages:
            A string specifying the species coverages. Unspecified species are
            set to zero. Can only be specified for interfaces.
        """
        self._t = temperature
        self._p = pressure
        self._rho = density
        self._x = mole_fractions
        self._y = mass_fractions
        self._c = coverages
        self._m = solute_molalities

    def build(self, ph):
        st = ph.addChild('state')
        if self._t: addFloat(st, 'temperature', self._t, defunits = 'K')
        if self._p: addFloat(st, 'pressure', self._p, defunits = _upres)
        if self._rho: addFloat(st, 'density', self._rho, defunits = _umass+'/'+_ulen+'3')
        if self._x: st.addChild('moleFractions', self._x)
        if self._y: st.addChild('massFractions', self._y)
        if self._c: st.addChild('coverages', self._c)
        if self._m: st.addChild('soluteMolalities', self._m)


class phase(object):
    """Base class for phases of matter."""

    def __init__(self,
                 name = '',
                 dim = 3,
                 elements = '',
                 species = '',
                 note = '',
                 reactions = 'none',
                 initial_state = None,
                 options = []):
        """
        :param name:
            A string to identify the phase. Must be unique among the phase
            names within the file.
        :param elements:
            The elements. A string of element symbols.
        :param species:
            The species. A string or sequence of strings in the format
            described in :ref:`sec-defining-species`.
        :param note:
            A user-defined comment. Not evaluated by Cantera itself.
        :param reactions:
            The homogeneous reactions. If omitted, no reactions will be
            included. A string or sequence of strings in the format described
            in :ref:`sec-declaring-reactions`. This field is not allowed for
            stoichiometric_solid and stoichiometric_liquid entries.
        :param kinetics:
            The kinetics model. Optional; if omitted, the default model for the
            phase type will be used.
        :param transport:
            The transport property model. Optional. If omitted, transport
            property calculation will be disabled.
        :param initial_state:
            Initial thermodynamic state, specified with an embedded state entry.
        :param options:
            Special processing options. Optional.
        """

        self._name = name
        self._dim = dim
        self._el = elements
        self._sp = []
        self._rx = []
        self._comment = note

        if isinstance(options, str):
            self._options = [options]
        else:
            self._options = options

        self.debug = 0
        if 'debug' in options:
            self.debug = 1

        #--------------------------------
        #        process species
        #--------------------------------

        # if a single string is entered, make it a list
        if isinstance(species, str):
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
            icolon = sp.find(':')
            if icolon > 0:
                #datasrc, spnames = sp.split(':')
                datasrc = sp[:icolon].strip()
                spnames = sp[icolon+1:]
                self._sp.append((datasrc+'.xml', spnames))
            else:
                spnames = sp
                self._sp.append(('', spnames))

            # strip the commas, and make the list of species names
            # 10/31/03: commented out the next line, so that species names may contain commas
            #sptoks = spnames.replace(',',' ').split()
            sptoks = spnames.split()

            for s in sptoks:
                # check for stray commas
                if s != ',':
                    if s[0] == ',': s = s[1:]
                    if s[-1] == ',': s = s[:-1]

                    if s != 'all' and s in self._spmap:
                        raise CTI_Error('Multiply-declared species '+s+' in phase '+self._name)
                    self._spmap[s] = self._dim

        self._rxns = reactions

        # check that species have been declared
        if len(self._spmap) == 0:
            raise CTI_Error('No species declared for phase '+self._name)

        # and that only one species is declared if it is a pure phase
        if self.is_pure() and len(self._spmap) > 1:
            raise CTI_Error('Stoichiometric phases must  declare exactly one species, \n'+
                               'but phase '+self._name+' declares '+str(len(self._spmap))+'.')

        self._initial = initial_state

        # add this phase to the global phase list
        global _phases
        _phases.append(self)


    def is_ideal_gas(self):
        """True if the entry represents an ideal gas."""
        return 0

    def is_pure(self):
        return 0

    def has_species(self, s):
        """Return 1 is a species with name 's' belongs to the phase,
        or 0 otherwise."""
        if s in self._spmap: return 1
        return 0

    def conc_dim(self):
        """Concentration dimensions. Used in computing the units for reaction
        rate coefficients."""
        return (1, -self._dim)


    def buildrxns(self, p):

        if isinstance(self._rxns, str):
            self._rxns = [self._rxns]

        # for each reaction string, check whether or not the reactions
        # are imported or defined locally. If imported, the string
        # contains a colon (:)
        for r in self._rxns:
            icolon = r.find(':')
            if icolon > 0:
                #datasrc, rnum = r.split(':')
                datasrc = r[:icolon].strip()
                rnum = r[icolon+1:]
                self._rx.append((datasrc+'.xml', rnum))
            else:
                rnum = r
                self._rx.append(('', rnum))

        for r in self._rx:
            datasrc = r[0]
            ra = p.addChild('reactionArray')
            ra['datasrc'] = datasrc+'#reaction_data'
            rk = None
            if 'skip_undeclared_species' in self._options:
                rk = ra.addChild('skip')
                rk['species'] = 'undeclared'

            if 'skip_undeclared_third_bodies' in self._options:
                if not rk:
                    rk = ra.addChild('skip')
                rk['third_bodies'] = 'undeclared'

            rtoks = r[1].split()
            if rtoks[0] != 'all':
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
        ph['dim'] = repr(self._dim)

        # ------- error tests -------
        #err = ph.addChild('validation')
        #err.addChild('duplicateReactions','halt')
        #err.addChild('thermo','warn')

        e = ph.addChild('elementArray',self._el)
        e['datasrc'] = 'elements.xml'
        for s in self._sp:
            datasrc, names = s
            sa = ph.addChild('speciesArray',names)
            sa['datasrc'] = datasrc+'#species_data'

            if 'skip_undeclared_elements' in self._options:
                sk = sa.addChild('skip')
                sk['element'] = 'undeclared'

        if self._rxns != 'none':
            self.buildrxns(ph)

        #self._eos.build(ph)
        if self._initial:
            self._initial.build(ph)

        if self._comment:
            ph.addChild('note',self._comment)

        thermo = ph.addChild('thermo')
        if 'allow_discontinuous_thermo' in self._options:
            thermo['allow_discontinuities'] = 'true'

        return ph


class ideal_gas(phase):
    """An ideal gas mixture."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 reactions = 'none',
                 kinetics = 'GasKinetics',
                 transport = 'None',
                 initial_state = None,
                 options = []):
        """
        The parameters correspond to those of :class:`.phase`, with the
        following modifications:

        :param kinetics:
            The kinetics model. Usually this field is omitted, in which case
            kinetics model GasKinetics, appropriate for reactions in ideal gas
            mixtures, is used.
        :param transport:
            The transport property model. One of the strings ``'none'``,
            ``'multi'``, or ``'mix'``. Default: ``'none'``.
        """

        phase.__init__(self, name, 3, elements, species, note, reactions,
                       initial_state, options)
        self._pure = 0
        self._kin = kinetics
        self._tr = transport
        if self.debug:
            print('Read ideal_gas entry '+self._name)
            try:
                print('in file '+__name__)
            except:
                pass



    def build(self, p):
        ph = phase.build(self, p)
        ph.child('thermo')['model'] = 'IdealGas'
        k = ph.addChild("kinetics")
        k['model'] = self._kin
        t = ph.addChild('transport')
        t['model'] = self._tr

    def is_ideal_gas(self):
        return 1


class stoichiometric_solid(phase):
    """
    A solid compound or pure element. Stoichiometric solid phases contain
    exactly one species, which always has unit activity. The solid is assumed
    to have constant density. Therefore the rates of reactions involving these
    phases do not contain any concentration terms for the (one) species in the
    phase, since the concentration is always the same."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 density = None,
                 transport = 'None',
                 initial_state = None,
                 options = []):
        """
        See :class:`.phase` for descriptions of the parameters.
        """

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._dens = density
        self._pure = 1
        if self._dens is None:
            raise CTI_Error('density must be specified.')
        self._tr = transport

    def conc_dim(self):
        """A stoichiometric solid always has unit activity, so the
        generalized concentration is 1 (dimensionless)."""
        return (0,0)

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child('thermo')
        e['model'] = 'StoichSubstance'
        addFloat(e, 'density', self._dens, defunits = _umass+'/'+_ulen+'3')
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'


class stoichiometric_liquid(stoichiometric_solid):
    """
    An incompressible stoichiometric liquid. Currently, there is no
    distinction between stoichiometric liquids and solids.
    """
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None,
                 options = []):
        """
        See :class:`.phase` for descriptions of the parameters.
        """
        stoichiometric_solid.__init__(self, name, elements,
                                      species, note, density, transport,
                                      initial_state, options)


class metal(phase):
    """A metal."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None,
                 options = []):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._dens = density
        self._pure = 0
        self._tr = transport

    def conc_dim(self):
        return (0,0)

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child("thermo")
        e['model'] = 'Metal'
        addFloat(e, 'density', self._dens, defunits = _umass+'/'+_ulen+'3')
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'

class semiconductor(phase):
    """A semiconductor."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 density = -1.0,
                 bandgap = 1.0 * eV,
                 effectiveMass_e = 1.0 * ElectronMass,
                 effectiveMass_h = 1.0 * ElectronMass,
                 transport = 'None',
                 initial_state = None,
                 options = []):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._dens = density
        self._pure = 0
        self._tr = transport
        self._emass = effectiveMass_e
        self._hmass = effectiveMass_h
        self._bandgap = bandgap

    def conc_dim(self):
        return (1,-3)

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child("thermo")
        e['model'] = 'Semiconductor'
        addFloat(e, 'density', self._dens, defunits = _umass+'/'+_ulen+'3')
        addFloat(e, 'effectiveMass_e', self._emass, defunits = _umass)
        addFloat(e, 'effectiveMass_h', self._hmass, defunits = _umass)
        addFloat(e, 'bandgap', self._bandgap, defunits = 'eV')
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'


class incompressible_solid(phase):
    """An incompressible solid."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None,
                 options = []):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._dens = density
        self._pure = 0
        if self._dens < 0.0:
            raise CTI_Error('density must be specified.')
        self._tr = transport

    def conc_dim(self):
        return (1,-3)

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child("thermo")
        e['model'] = 'Incompressible'
        addFloat(e, 'density', self._dens, defunits = _umass+'/'+_ulen+'3')
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'


class lattice(phase):
    def __init__(self, 
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 reactions = 'none',
                 transport = 'None',
                 initial_state = None,
                 options = [],
                 site_density = -1.0,
                 vacancy_species = ''):
        phase.__init__(self, name, 3, elements, species, note, 'none',
                        initial_state, options)
        self._tr = transport
        self._n = site_density
        self._vac = vacancy_species
        self._species = species
        if name == '':
            raise CTI_Error('sublattice name must be specified')
        if species == '':
            raise CTI_Error('sublattice species must be specified')
        if site_density < 0.0:
            raise CTI_Error('sublattice '+name
                            +' site density must be specified')

    def build(self,p, visible = 0):
        #if visible == 0:
        #    return
        ph = phase.build(self, p)
        e = ph.child('thermo')
        e['model'] = 'Lattice'
        addFloat(e, 'site_density', self._n, defunits = _umol+'/'+_ulen+'3')
        if self._vac:
            e.addChild('vacancy_species',self._vac)
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'

class lattice_solid(phase):
    """A solid crystal consisting of one or more sublattices."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 lattices = [],
                 transport = 'None',
                 initial_state = None,
                 options = []):

        # find elements
        elist = []
        for lat in lattices:
            e = lat._el.split()
            for el in e:
                if not el in elist:
                    elist.append(el)
        elements = ' '.join(elist)

        # find species
        slist = []
        for lat in lattices:
            _sp = ""
            for spp in lat._species:
                _sp = _sp + spp
            s = _sp.split()
            for sp in s:
                if not sp in slist:
                    slist.append(sp)
        species = ' '.join(slist)

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._lattices = lattices
        if lattices == []:
            raise CTI_Error('One or more sublattices must be specified.')
        self._pure = 0
        self._tr = transport

    def conc_dim(self):
        return (0,0)

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child("thermo")
        e['model'] = 'LatticeSolid'

        if self._lattices:
            lat = e.addChild('LatticeArray')
            for n in self._lattices:
                n.build(lat, visible = 1)

        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'



class liquid_vapor(phase):
    """A fluid with a complete liquid/vapor equation of state.
    This entry type selects one of a set of predefined fluids with
    built-in liquid/vapor equations of state. The substance_flag
    parameter selects the fluid. See purefluids.py for the usage
    of this entry type."""

    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 substance_flag = 0,
                 initial_state = None,
                 options = []):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._subflag = substance_flag
        self._pure = 1


    def conc_dim(self):
        return (0,0)

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child("thermo")
        e['model'] = 'PureFluid'
        e['fluid_type'] = repr(self._subflag)
        k = ph.addChild("kinetics")
        k['model'] = 'none'

class RedlichKwongMFTP(phase):
    """A multi-component fluid model for non-ideal gas fluids.
        """

    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 initial_state = None,
                 activity_coefficients = None,
                 transport = 'None',
                 options = []):

        phase.__init__(self,name, 3, elements, species, note, 'none',
                       initial_state,options)
        self._pure = 0
        self._tr = transport
        self._activityCoefficients = activity_coefficients

    def conc_dim(self):
        return (0,0)

    def build(self, p):
        ph = phase.build(self,p)
        e = ph.child("thermo")
        e['model'] = 'RedlichKwongMFTP'
        if self._activityCoefficients:
            a = e.addChild("activityCoefficients")
            if isinstance(self._activityCoefficients, activityCoefficients):
                self._activityCoefficients.build(a)
            else:
                na = len(self._activityCoefficients)
                for n in range(na):
                    self._activityCoefficients[n].build(a)

        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'


class redlich_kwong(phase):
    """A fluid with a complete liquid/vapor equation of state.
    This entry type selects one of a set of predefined fluids with
    built-in liquid/vapor equations of state. The substance_flag
    parameter selects the fluid. See purefluids.py for the usage
    of this entry type."""

    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 substance_flag = 7,
                 initial_state = None,
                 Tcrit = 1.0,
                 Pcrit = 1.0,
                 options = []):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._subflag = 7
        self._pure = 1
        self._tc = 1
        self._pc = 1

    def conc_dim(self):
        return (0,0)

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child("thermo")
        e['model'] = 'PureFluid'
        e['fluid_type'] = repr(self._subflag)
        addFloat(e, 'Tc', self._tc, defunits = "K")
        addFloat(e, 'Pc', self._pc, defunits = "Pa")
        addFloat(e, 'MolWt', self._mw, defunits = _umass+"/"+_umol)
        ph.addChild("kinetics")
        k['model'] = 'none'



class ideal_interface(phase):
    """A chemically-reacting ideal surface solution of multiple species."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 reactions = 'none',
                 site_density = 0.0,
                 phases = [],
                 kinetics = 'Interface',
                 transport = 'None',
                 initial_state = None,
                 options = []):
        """
        The parameters correspond to those of :class:`.phase`, with the
        following modifications:

        :param reactions:
            The heterogeneous reactions at this interface. If omitted, no
            reactions will be included. A string or sequence of strings in the
            format described in :ref:`sec-declaring-reactions`.
        :param site_density:
            The number of adsorption sites per unit area.
        :param phases:
            A string listing the bulk phases that participate in reactions
            at this interface.
        """
        self._type = 'surface'
        phase.__init__(self, name, 2, elements, species, note, reactions,
                       initial_state, options)
        self._pure = 0
        self._kin = kinetics
        self._tr = transport
        self._phases = phases
        self._sitedens = site_density

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child("thermo")
        e['model'] = 'Surface'
        addFloat(e, 'site_density', self._sitedens, defunits = _umol+'/'+_ulen+'2')
        k = ph.addChild("kinetics")
        k['model'] = self._kin
        t = ph.addChild('transport')
        t['model'] = self._tr
        p = ph.addChild('phaseArray',self._phases)


    def conc_dim(self):
        return (1, -2)


class edge(phase):
    """A 1D boundary between two surface phases."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 reactions = 'none',
                 site_density = 0.0,
                 phases = [],
                 kinetics = 'Edge',
                 transport = 'None',
                 initial_state = None,
                 options = []):

        self._type = 'edge'
        phase.__init__(self, name, 1, elements, species, note, reactions,
                       initial_state, options)
        self._pure = 0
        self._kin = kinetics
        self._tr = transport
        self._phases = phases
        self._sitedens = site_density

    def build(self, p):
        ph = phase.build(self, p)
        e = ph.child("thermo")
        e['model'] = 'Edge'
        addFloat(e, 'site_density', self._sitedens, defunits = _umol+'/'+_ulen)
        k = ph.addChild("kinetics")
        k['model'] = self._kin
        t = ph.addChild('transport')
        t['model'] = self._tr
        p = ph.addChild('phaseArray',self._phases)


    def conc_dim(self):
        return (1, -1)


## class binary_salt_parameters:
##     def __init__(self,
##                  cation = "",
##                  anion = "",
##                  beta0 = None,
##                  beta1 = None,
##                  beta2 = None,
##                  Cphi = None,
##                  Alpha1 = -1.0):
##         self._cation = cation
##         self._anion = anion
##         self._beta0 = beta0
##         self._beta1 = beta1
##         self._Cphi = Cphi
##         self._Alpha1 = Alpha1

##     def build(self, a):
##         s = a.addChild("binarySaltParameters")
##         s["cation"] = self._cation
##         s["anion"] = self._anion
##         s.addChild("beta0", self._beta0)
##         s.addChild("beta1", self._beta1)
##         s.addChild("beta2", self._beta2)
##         s.addChild("Cphi", self._Cphi)
##         s.addChild("Alpha1", self._Alpha1)

## class theta_anion:
##     def __init__(self,
##                  anions = None,
##                  theta = 0.0):
##         self._anions = anions
##         self._theta = theta

##     def build(self, a):
##         s = a.addChild("thetaAnion")
##         s["anion1"] = self._anions[0]
##         s["anion2"] = self._anions[1]
##         s.addChild("Theta", self._theta)

## class psi_common_cation:
##     def __init__(self,
##                  anions = None,
##                  cation = '',
##                  theta = 0.0,
##                  psi = 0.0):
##         self._anions = anions
##         self._cation = cation
##         self._theta = theta
##         self._psi = psi

##     def build(self, a):
##         s = a.addChild("psiCommonCation")
##         s["anion1"] = self._anions[0]
##         s["anion2"] = self._anions[1]
##         s["cation"] = self._cation
##         s.addChild("Theta", self._theta)
##         s.addChild("Psi", self._psi)

## class psi_common_anion:
##     def __init__(self,
##                  anion = '',
##                  cations = None,
##                  theta = 0.0,
##                  psi = 0.0):
##         self._anion = anion
##         self._cations = cations
##         self._theta = theta
##         self._psi = psi

##     def build(self, a):
##         s = a.addChild("psiCommonAnion")
##         s["anion1"] = self._cations[0]
##         s["anion2"] = self._cations[1]
##         s["cation"] = self._anion
##         s.addChild("Theta", self._theta)
##         s.addChild("Psi", self._psi)


## class theta_cation:
##     def __init__(self,
##                  cations = None,
##                  theta = 0.0):
##         self._cations = cations
##         self._theta = theta

##     def build(self, a):
##         s = a.addChild("thetaCation")
##         s["cation1"] = self._anions[0]
##         s["cation2"] = self._anions[1]
##         s.addChild("Theta", self._theta)

## class pitzer:
##     def __init__(self,
##                  temp_model = "",
##                  A_Debye = "",
##                  default_ionic_radius = -1.0,

## class electrolyte(phase):
##     """An electrolye solution obeying the HMW model."""
##     def __init__(self,
##                  name = '',
##                  elements = '',
##                  species = '',
##                  note = '',
##                  transport = 'None',
##                  initial_state = None,
##                  solvent = '',
##                  standard_concentration = '',
##                  activity_coefficients = None,
##                  options = []):

##         phase.__init__(self, name, 3, elements, species, note, 'none',
##                        initial_state, options)
##         self._pure = 0
##         self._solvent = solvent
##         self._stdconc = standard_concentration

##     def conc_dim(self):
##         return (1,-3)

##     def build(self, p):
##         ph = phase.build(self, p)
##         e = ph.child("thermo")
##         sc = e.addChild("standardConc")
##         sc['model'] = self._stdconc
##         e['model'] = 'HMW'
##         e.addChild("activity_coefficients")

##         addFloat(e, 'density', self._dens, defunits = _umass+'/'+_ulen+'3')
##         if self._tr:
##             t = ph.addChild('transport')
##             t['model'] = self._tr
##         k = ph.addChild("kinetics")
##         k['model'] = 'none'


#-------------------------------------------------------------------

# falloff parameterizations

class Troe(object):
    """The Troe falloff function."""
    def __init__(self, A = 0.0, T3 = 0.0, T1 = 0.0, T2 = -999.9):
        """
        Parameters: *A*, *T3*, *T1*, *T2*. These must be entered as pure
        numbers with no attached dimensions.
        """
        if T2 != -999.9:
            self._c = (A, T3, T1, T2)
        else:
            self._c = (A, T3, T1)

    def build(self, p):
        s = ''
        for num in self._c:
            s += '%g ' % num
        f = p.addChild('falloff', s)
        f['type'] = 'Troe'


class SRI(object):
    """ The SRI falloff function."""
    def __init__(self, A = 0.0, B = 0.0, C = 0.0, D = -999.9, E=-999.9):
        """
        Parameters: *A*, *B*, *C*, *D*, *E*. These must be entered as
        pure numbers without attached dimensions.
        """
        if D != -999.9 and E != -999.9:
            self._c = (A, B, C, D, E)
        else:
            self._c = (A, B, C)

    def build(self, p):
        s = ''
        for num in self._c:
            s += '%g ' % num
        f = p.addChild('falloff', s)
        f['type'] = 'SRI'


class Lindemann(object):
    """The Lindemann falloff function."""
    def __init__(self):
        """ This falloff function takes no parameters."""
        pass
    def build(self, p):
        f = p.addChild('falloff')
        f['type'] = 'Lindemann'


#get_atomic_wts()
validate()

def convert(filename, outName=None):
    import os, sys
    base = os.path.basename(filename)
    root, _ = os.path.splitext(base)
    dataset(root)
    try:
        with open(filename, 'rU') as f:
            code = compile(f.read(), filename, 'exec')
            exec(code)
    except SyntaxError as err:
        # Show more context than the default SyntaxError message
        # to help see problems in multi-line statements
        text = open(filename, 'rU').readlines()
        print('%s in "%s" on line %i:\n' % (err.__class__.__name__,
                                            err.filename,
                                            err.lineno))
        print('|  Line |')
        for i in range(max(err.lineno-6, 0),
                       min(err.lineno+3, len(text))):
            print('| % 5i |' % (i+1), text[i].rstrip())
            if i == err.lineno-1:
                print(' '* (err.offset+9) + '^')
        print()
        sys.exit(3)
    except TypeError as err:
        import traceback

        text = open(filename, 'rU').readlines()
        tb = traceback.extract_tb(sys.exc_info()[2])
        lineno = tb[-1][1]

        print('%s on line %i of %s:' % (err.__class__.__name__, lineno, filename))
        print(err)
        print('\n| Line |')

        for i in range(max(lineno-6, 0),
                       min(lineno+3, len(text))):
            if i == lineno-1:
                print('> % 4i >' % (i+1), text[i].rstrip())
            else:
                print('| % 4i |' % (i+1), text[i].rstrip())

        sys.exit(4)

    write(outName)

if __name__ == "__main__":
    import sys
    if len(sys.argv) not in (2,3):
        raise ValueError('Incorrect number of command line arguments.')
    convert(*sys.argv[1:])
