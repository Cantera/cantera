#!/usr/bin/env python

# This file is part of Cantera. See License.txt in the top-level directory or
# at htts://www.cantera.org/license.txt for license and copyright information.

"""
cti2yaml.py: Convert legacy CTI input files to YAML

Usage:
    python ck2yaml.py mech.cti

This will produce the output file 'mech.yaml'.
"""

from __future__ import print_function

import sys
import math
import numpy as np

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

# Python 2/3 compatibility
try:
  basestring
except NameError:
  basestring = str


def _printerr(*args):
    # All debug and error output should go to stderr
    print(*args, file=sys.stderr)


class InputError(Exception):
    """
    Exception raised if an error is encountered while parsing the input file.
    """
    def __init__(self, msg):
        _printerr('\n\n***** Error parsing input file *****\n\n')
        _printerr(msg)
        _printerr()

BlockMap = yaml.comments.CommentedMap

def FlowMap(*args, **kwargs):
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m

def FlowList(*args, **kwargs):
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst

def represent_float(self, data):
    # type: (Any) -> Any
    if data != data:
        value = u'.nan'
    elif data == self.inf_value:
        value = u'.inf'
    elif data == -self.inf_value:
        value = u'-.inf'
    else:
        if data == 0:
            value = u'0.0'
        elif 0.01 <= abs(data) < 10000:
            value = np.format_float_positional(data, trim='0')
        else:
            value = np.format_float_scientific(data, trim='0')

    return self.represent_scalar(u'tag:yaml.org,2002:float', value)

yaml.RoundTripRepresenter.add_representer(float, represent_float)


def applyUnits(value):
    if isinstance(value, (float, int)):
        return value
    else:
        if value[0] == 0:
            return '0.0 {}'.format(value[1])
        elif 0.01 <= abs(value[0]) < 10000:
            return '{} {}'.format(
                np.format_float_positional(value[0], trim='0'), value[1])
        else:
            return '{} {}'.format(
                np.format_float_scientific(value[0], trim='0'), value[1])


# map of old CTI/XML names with the new YAML names
_newNames = {
    'GasKinetics': 'gas',
    'Mix': 'mixture-averaged',
    'Multi': 'multicomponent'
}

# constants that can be used in .cti files
OneAtm = 1.01325e5
OneBar = 1.0e5
# Conversion from eV to J/kmol (electronCharge * Navrog)
eV = 9.64853364595687e7
# Electron Mass in kg
ElectronMass = 9.10938291e-31

# default units
_ulen = 'm'
_umol = 'kmol'
_umass = 'kg'
_utime = 's'
_ue = 'J/kmol'
_uenergy = 'J'
_upres = 'Pa'

# default std state pressure
_pref = OneAtm

_name = 'noname'

# these lists store top-level entries
_elements = []
_species = []
_speciesnames = []
_phases = []
_reactions = []

# default for Motz & Wise correction
_motz_wise = None

def enable_motz_wise():
    global _motz_wise
    _motz_wise = True

def disable_motz_wise():
    global _motz_wise
    _motz_wise = False

def validate(species = 'yes', reactions = 'yes'):
    pass

def dataset(nm):
    "Set the dataset name. Invoke this to change the name of the XML file."
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


def addFloat(x, nm, val, fmt='', defunits=''):
    """
    Add a child element to XML element x representing a
    floating-point number.
    """
    u = ''
    s = ''
    if isinstance(val, (float, int)):
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
        try:
            d[b[0]] = int(b[1])
        except ValueError:
            d[b[0]] = float(b[1])
    return d


class element(object):
    """ An atomic element or isotope. """
    def __init__(self, symbol = '',
                 atomic_mass = 0.01,
                 atomic_number = None):
        """
        :param symbol:
            The symbol for the element or isotope.
        :param atomic_mass:
            The atomic mass in amu.
        """
        self._sym = symbol
        self._atw = atomic_mass
        self._num = atomic_number
        _elements.append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('symbol', node._sym),
                        ('atomic-weight', node._atw)])
        if node._num is not None:
            out['atomic-number'] = node._num
        return representer.represent_dict(out)


class species(object):
    """A constituent of a phase or interface."""

    def __init__(self,
                 name = 'missing name!',
                 atoms = '',
                 note = '',
                 thermo = None,
                 transport = None,
                 charge = None,
                 size = 1.0):
        """
        :param name:
            The species name (or formula). The name may be arbitrarily long,
            although usually a relatively short, abbreviated name is most
            convenient. Required parameter.
        :param atoms:
            The atomic composition, specified by a string containing
            space-delimited ``<element>:<atoms>`` pairs. The number of atoms may be
            either an integer or a floating-point number.
        :param note:
            A user-defined comment. Not evaluated by Cantera itself.
        :param thermo:
            The parameterization to use to compute the reference-state
            thermodynamic properties. This must be one of the entry types
            described in `Thermodynamic Property Models
            <https://cantera.org/science/science-species.html#sec-thermo-models>`__.
            To specify multiple parameterizations, each for a different temperature range,
            group them in parentheses.
        :param transport:
            An entry specifying parameters to compute this species'
            contribution to the transport properties. This must be one of the
            entry types described in `Species Transport Coefficients
            <https://cantera.org/science/science-species.html#species-transport-coefficients>`__,
            and must be consistent with the transport model of the phase into which
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
        if charge is not None and 'E' not in self._atoms:
            self._atoms['E'] = -charge
        self._size = size
        self._comment = note

        if isinstance(thermo, (list, tuple)):
            if isinstance(thermo[0], (NASA, NASA9)):
                self._thermo = MultiNASA(thermo)
        elif isinstance(thermo[0], (NASA, NASA9)):
            self._thermo = MultiNASA([thermo])
        elif thermo is not None:
            self._thermo = thermo
        else:
            self._thermo = const_cp()

        self._transport = transport

        _species.append(self)
        _speciesnames.append(name)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('name', node._name),
                        ('composition', FlowMap(node._atoms.items()))])
        if node._size != 1:
            out['sites'] = node._size
        out['thermo'] = node._thermo
        if node._transport:
            out['transport'] = node._transport
        if node._comment:
            out['note'] = node._comment
        return representer.represent_dict(out)


class thermo(object):
    """Base class for species standard-state thermodynamic properties."""
    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)


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
    def __init__(self, Trange = (0.0, 0.0), coeffs = [], p0 = None):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required.
        :param coeffs:
            List of seven coefficients :math:`(a_0, \ldots , a_6)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value is used, which is set by the ``standard_pressure``
            directive.
        """
        self._t = Trange
        self._pref = p0
        if len(coeffs) != 7:
            raise InputError('NASA coefficient list must have length = 7')
        self._coeffs = coeffs


class NASA9(thermo):
    """NASA9 polynomial parameterization for a single temperature region."""

    def __init__(self, Trange = (0.0, 0.0),
                 coeffs = [], p0 = -1.0):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required.
        :param coeffs:
            List of nine coefficients :math:`(a_0, \ldots , a_8)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value is used, which is set by the ``standard_pressure``
            directive.
        """
        self._t = Trange         # Range of the polynomial representation
        self._pref = p0          # Reference pressure
        if len(coeffs) != 9:
            raise InputError('NASA9 coefficient list must have length = 9')
        self._coeffs = coeffs


class MultiNASA(thermo):
    def __init__(self, regions):
        regions = sorted(regions, key=lambda r: r._t[0])
        self.pref = regions[0]._pref or _pref
        self.Tranges = [regions[0]._t[0]]
        self.model = 'NASA7' if len(regions[0]._coeffs) == 7 else 'NASA9'
        self.data = []
        for r in regions:
            self.Tranges.append(r._t[1])
            self.data.append(r._coeffs)

    def get_yaml(self, out):
        out['model'] = self.model
        if self.pref != OneAtm:
            out['reference-pressure'] = self.pref

        out['temperature-ranges'] = FlowList(self.Tranges)
        out['data'] = [FlowList(coeffs) for coeffs in self.data]


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
            raise InputError('Shomate coefficient list must have length = 7')
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
        c = t.addChild('const_cp')
        if self._t[0] >= 0.0: c['Tmin'] = repr(self._t[0])
        if self._t[1] >= 0.0: c['Tmax'] = repr(self._t[1])
        energy_units = _uenergy+'/'+_umol
        addFloat(c,'t0',self._c[0], defunits = 'K')
        addFloat(c,'h0',self._c[1], defunits = energy_units)
        addFloat(c,'s0',self._c[2], defunits = energy_units+'/K')
        addFloat(c,'cp0',self._c[3], defunits = energy_units+'/K')


class gas_transport(object):
    """
    Species-specific Transport coefficients for gas-phase transport models.
    """
    def __init__(self, geom,
                 diam = 0.0, well_depth = 0.0, dipole = 0.0,
                 polar = 0.0, rot_relax = 0.0, acentric_factor = None,
                 disp_coeff = 0.0, quad_polar = 0.0):
        """
        :param geom:
            A string specifying the molecular geometry. One of ``atom``,
            ``linear``, or ``nonlinear``. Required.
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
        :param w_ac:
            Pitzer's acentric factor.  Dimensionless.
            Default: 0.0
        :param disp_coeff:
            The dispersion coefficient in A^5
            Default: 0.0
        :param quad_polar:
            The quadrupole polarizability
            Default: 0.0
        """
        self._geom = geom
        self._diam = diam
        self._well_depth = well_depth
        self._dipole = dipole
        self._polar = polar
        self._rot_relax = rot_relax
        self._w_ac = acentric_factor
        self._disp_coeff = disp_coeff
        self._quad_polar = quad_polar

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('model', 'gas'),
                        ('geometry', node._geom),
                        ('diameter', node._diam),
                        ('well-depth', node._well_depth)])
        if node._dipole:
            out['dipole'] = node._dipole
        if node._polar:
            out['polarizability'] = node._polar
        if node._rot_relax:
            out['rotational-relaxation'] = node._rot_relax
        if node._w_ac:
            out['acentric-factor'] = node._w_ac
        if node._disp_coeff:
            out['dispersion-coefficient'] = node._disp_coeff
        if node._quad_polar:
            out['quadrupole-polarizability'] = node._quad_polar
        return representer.represent_dict(out)


class rate_expression(object):
    pass

class Arrhenius(rate_expression):
    def __init__(self,
                 A = 0.0,
                 b = 0.0,
                 E = 0.0,
                 coverage = []):
        """
        :param A:
            The pre-exponential coefficient. Required input. If entered without
            units, the units will be computed considering all factors that
            affect the units. The resulting units string is written to the CTML
            file individually for each reaction pre-exponential coefficient.
        :param b:
            The temperature exponent. Dimensionless. Default: 0.0.
        :param E:
            Activation energy. Default: 0.0.
        :param coverage: For a single coverage dependency, a list with four
            elements: the species name followed by the three coverage
            parameters. For multiple coverage dependencies, a list of lists
            containing the individual sets of coverage parameters. Only used for
            surface and edge reactions.
        """

        self._c = [A, b, E]

        if coverage:
            if isinstance(coverage[0], basestring):
                self._cov = [coverage]
            else:
                self._cov = coverage
            for cov in self._cov:
                if len(cov) != 4:
                    raise InputError("Incorrect number of coverage parameters")
        else:
            self._cov = None

    @classmethod
    def to_yaml(cls, representer, node):
        out = FlowMap([('A', applyUnits(node._c[0])),
                       ('b', applyUnits(node._c[1])),
                       ('Ea', applyUnits(node._c[2]))])
        return representer.represent_dict(out)


class stick(Arrhenius):
    def __init__(self, *args, **kwargs):
        """
        :param motz_wise: 'True' if the Motz & Wise correction should be used,
            'False' if not. If unspecified, use the mechanism default (set using
            the functions `enable_motz_wise` or `disable_motz_wise`).
        """
        self.motz_wise = kwargs.pop('motz_wise', None)
        Arrhenius.__init__(self, *args, **kwargs)

    def build(self, p, name=''):
        a = p.addChild('Arrhenius')
        a['type'] = 'stick'
        ngas = len(self.gas_species)
        if ngas != 1:
            raise InputError("Sticking probabilities can only be used for "
                "reactions with one gas-phase reactant, but this reaction has "
                + str(ngas) + ': ' + str(self.gas_species))

        a['species'] = self.gas_species[0]
        if self.motz_wise is not None:
            a['motz_wise'] = str(self.motz_wise).lower()
        Arrhenius.build(self, p, name, a)

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
        r"""
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function :math:`A T^b exp(-E/\hat{R}T)`.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param order:
            Override the default reaction orders implied by the reactant
            stoichiometric coefficients. Given as a string of key:value pairs,
            e.g., ``"CH4:0.25 O2:1.5"``.
        :param options: Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
            May be one or more (as a list) of the
            following: ``'duplicate'``, ``'negative_A'``,`` 'negative_orders'``,
            ``'nonreactant_orders'``.
        """
        self._e = equation
        self._order = getPairs(order)
        self._num = len(_reactions)+1

        if isinstance(options, str):
            self._options = [options]
        else:
            self._options = options

        if isinstance(kf, (list, tuple)):
            kf = Arrhenius(*kf)
        self._kf = kf
        self._type = 'elementary'
        _reactions.append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)

    def get_yaml(self, out):
        out['equation'] = self._e
        out.yaml_add_eol_comment('Reaction {}'.format(self._num), 'equation')
        if self._type != 'elementary':
            out['type'] = self._type

        if self._type in ('elementary', 'three-body', 'edge', 'surface'):
            out['rate'] = self._kf

        if 'duplicate' in self._options:
            out['duplicate'] = True
        if 'negative_A' in self._options:
            out['negative-A'] = True

        if self._order:
            out['orders'] = FlowMap(self._order.items())
        if 'negative_orders' in self._options:
            out['negative-orders'] = True
        if 'nonreactant_orders' in self._options:
            out['nonreactant-orders'] = True

    def build(self, p):
        # The default rate coefficient type is Arrhenius. If the rate
        # coefficient has been specified as a sequence of three
        # numbers, then create a new Arrhenius instance for it;
        # otherwise, just use the supplied instance.
        nm = ''

        if self._type == 'edge' or self._type == 'surface':
            if self._beta > 0:
                electro = kfnode.addChild('electrochem')
                electro['beta'] = repr(self._beta)


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
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies.
            The efficiencies for unspecified species are set to 1.0.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options: Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
        """
        reaction.__init__(self, equation, kf, id, '', options)
        self._type = 'three-body'
        self._eff = getPairs(efficiencies)

    def get_yaml(self, out):
        super(three_body_reaction, self).get_yaml(out)
        if self._eff:
            out['efficiencies'] = FlowMap(self._eff)


class pdep_reaction(reaction):
    """ Base class for falloff_reaction and chemically_activated_reaction """
    def __init__(self, equation, klow, khigh, efficiencies, falloff, id, options):
        super(pdep_reaction, self).__init__(equation, None, id, '', options)
        self._klow = klow
        self._khigh = khigh
        self._falloff = falloff
        self._eff = getPairs(efficiencies)

    def get_yaml(self, out):
        super(pdep_reaction, self).get_yaml(out)

        if isinstance(self._klow, (list, tuple)):
            self._klow = Arrhenius(*self._klow)
        if isinstance(self._khigh, (list, tuple)):
            self._khigh = Arrhenius(*self._khigh)

        out['low-rate'] = self._klow
        out['high-rate'] = self._khigh

        if self._falloff:
            self._falloff.get_yaml(out)

        if self._eff:
            out['efficiencies'] = FlowMap(self._eff)


class falloff_reaction(pdep_reaction):
    """ A gas-phase falloff reaction. """
    def __init__(self, equation, kf0, kf,
                 efficiencies='', falloff=None, id='', options=[]):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction in the high-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param kf0:
            The rate coefficient for the forward direction in the low-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
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
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
        """
        super(falloff_reaction, self).__init__(equation, kf0, kf, efficiencies,
            falloff, id, options)
        self._type = 'falloff'


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
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param kHigh:
            The rate coefficient for the forward direction in the high-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
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
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
        """
        super(chemically_activated_reaction, self).__init__(equation, kLow, 
            kHigh, efficiencies, falloff, id, options)
        self._type = 'chemically-activated'


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
        for p, A, b, Ea in args:
            self.pressures.append(p)
            self.arrhenius.append((A, b, Ea))

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

    def build(self, p):
        r = reaction.build(self, p)
        kfnode = r.child('rateCoeff')
        addFloat(kfnode, 'Tmin', self.Tmin)
        addFloat(kfnode, 'Tmax', self.Tmax)
        addFloat(kfnode, 'Pmin', self.Pmin)
        addFloat(kfnode, 'Pmax', self.Pmax)

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
    def __init__(self, equation='', kf=None, id='', order='', beta = 0.0,
                 options=[]):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function.
        :param sticking_prob:
            The reactive sticking probability for the forward direction. This
            can only be specified if there is only one bulk-phase reactant and
            it belongs to an ideal gas phase. If a sequence of three numbers is
            given, these will be interpreted as [A, b, E] in the modified
            Arrhenius function.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options:
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
        :param beta:
            Charge transfer coefficient: A number between 0 and 1 which, for a
            charge transfer reaction, determines how much of the electric
            potential difference between two phases is applied to the
            activiation energy of the fwd reaction.  The remainder is applied to
            the reverse reaction.
        """
        reaction.__init__(self, equation, kf, id, order, options)
        self._type = 'surface'
        self._beta = beta


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
            described in `Defining the Species
            <https://cantera.org/tutorials/cti/phases.html#defining-the-species>`__.
        :param note:
            A user-defined comment. Not evaluated by Cantera itself.
        :param reactions:
            The homogeneous reactions. If omitted, no reactions will be
            included. A string or sequence of strings in the format described
            in `Declaring the Reactions
            <https://cantera.org/tutorials/cti/phases.html#declaring-the-reactions>`__.
            This field is not allowed for ``stoichiometric_solid`` and
            ``stoichiometric_liquid`` entries.
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
        self._el = elements
        self._sp = []
        self._rx = []
        self._kin = None
        self._tr = None
        self._comment = note

        if isinstance(options, str):
            self._options = [options]
        else:
            self._options = options

        #--------------------------------
        #        process species
        #--------------------------------

        # if a single string is entered, make it a list
        if isinstance(species, str):
            species = [species]

        # for each species string, check whether or not the species
        # are imported or defined locally. If imported, the string
        # contains a colon (:)
        for sp in species:
            foundColon = False
            allLocal = True
            for token in sp.split():
                if ':' in sp:
                    foundColon = True
                if token not in _speciesnames:
                    allLocal = False

            if foundColon and not allLocal:
                icolon = sp.find(':')
                datasrc = sp[:icolon].strip()
                spnames = sp[icolon+1:]
                self._sp.append((datasrc + '.yaml', spnames.split()))
            else:
                spnames = sp
                self._sp.append(('species', spnames.split()))

        if isinstance(reactions, str):
            reactions = [reactions]

        # for each reaction string, check whether or not the reactions
        # are imported or defined locally. If imported, the string
        # contains a colon (:)
        for r in reactions:
            icolon = r.find(':')
            if icolon > 0:
                datasrc = r[:icolon].strip() + '.yaml'
                rnum = r[icolon+1:]
            else:
                datasrc = 'reactions'
                rnum = r
            if rnum == 'all' and 'skip_undeclared_species' in self._options:
                rnum = 'declared-species'
            self._rx.append((datasrc, rnum))
            if rnum not in ('all', 'declared-species'):
                _printerr("WARNING: Reactions specification '{}' not supported".format(
                    rnum))

        self._initial = initial_state

        # add this phase to the global phase list
        _phases.append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)

    def get_yaml(self, out):
        out['name'] = self._name
        out['thermo'] = self.thermo_model
        out['elements'] = FlowList(self._el.split())

        if len(self._sp) == 1 and self._sp[0][0] == 'species':
            # all local species
            out['species'] = FlowList(self._sp[0][1])
        else:
            out['species'] = [BlockMap([(sp[0], FlowList(sp[1]))])
                              for sp in self._sp]

        if 'skip_undeclared_elements' in self._options:
            out['skip-undeclared-elements'] = True

        if self._kin:
            out['kinetics'] = _newNames[self._kin]
            if len(self._rx) == 1 and self._rx[0][0] == 'reactions':
                out['reactions'] = self._rx[0][1]
            else:
                out['reactions'] = [BlockMap([(r[0], r[1])]) for r in self._rx]

        if self._tr:
            out['transport'] = _newNames[self._tr]

        if self._comment:
            out['note'] = self._comment

    def build(self, p):
        if self._initial:
            self._initial.build(ph)


class ideal_gas(phase):
    """An ideal gas mixture."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 note = '',
                 reactions = 'none',
                 kinetics = 'GasKinetics',
                 transport = None,
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
        self._kin = kinetics
        self._tr = transport
        self.thermo_model = 'ideal-gas'


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
        if self._dens is None:
            raise InputError('density must be specified.')
        self._tr = transport

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
        self._tr = transport

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
        self._tr = transport
        self._emass = effectiveMass_e
        self._hmass = effectiveMass_h
        self._bandgap = bandgap

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
                 density = None,
                 transport = 'None',
                 initial_state = None,
                 options = []):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._dens = density
        if self._dens is None:
            raise InputError('density must be specified.')
        self._tr = transport

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
                 site_density = None):
        phase.__init__(self, name, 3, elements, species, note, 'none',
                        initial_state, options)
        self._tr = transport
        self._n = site_density
        self._species = species
        if name == '':
            raise InputError('sublattice name must be specified')
        if species == '':
            raise InputError('sublattice species must be specified')
        if site_density is None:
            raise InputError('sublattice '+name
                            +' site density must be specified')

    def build(self,p, visible = 0):
        ph = phase.build(self, p)
        e = ph.child('thermo')
        e['model'] = 'Lattice'
        addFloat(e, 'site_density', self._n, defunits = _umol+'/'+_ulen+'3')
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
                _sp += spp
            s = _sp.split()
            for sp in s:
                if not sp in slist:
                    slist.append(sp)
        species = ' '.join(slist)

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self._lattices = lattices
        if lattices == []:
            raise InputError('One or more sublattices must be specified.')
        self._tr = transport

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
    parameter selects the fluid. See liquidvapor.cti and liquidvapor.py
    for the usage of this entry type."""

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
                 reactions = 'none',
                 kinetics = 'GasKinetics',
                 initial_state = None,
                 activity_coefficients = None,
                 transport = 'None',
                 options = []):

        phase.__init__(self,name, 3, elements, species, note, reactions,
                       initial_state,options)
        self._kin = kinetics
        self._tr = transport
        self._activityCoefficients = activity_coefficients

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
        if self._kin:
            k = ph.addChild("kinetics")
            k['model'] = self._kin


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
            format described in `Declaring the Reactions
            <https://cantera.org/tutorials/cti/phases.html#declaring-the-reactions>`__.
        :param site_density:
            The number of adsorption sites per unit area.
        :param phases:
            A string listing the bulk phases that participate in reactions
            at this interface.
        """
        self._type = 'surface'
        phase.__init__(self, name, 2, elements, species, note, reactions,
                       initial_state, options)
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

#-------------------------------------------------------------------

# falloff parameterizations

class Troe(object):
    """The Troe falloff function."""
    def __init__(self, A=0.0, T3=0.0, T1=0.0, T2=None):
        """
        Parameters: *A*, *T3*, *T1*, *T2*. These must be entered as pure
        numbers with no attached dimensions.
        """
        self.A = A
        self.T3 = T3
        self.T1 = T1
        self.T2 = T2

    def get_yaml(self, out):
        troe = FlowMap([('A', self.A), ('T3', self.T3), ('T1', self.T1)])
        if self.T2 is not None:
            troe['T2'] = self.T2
        out['Troe'] = troe


class SRI(object):
    """ The SRI falloff function."""
    def __init__(self, A=0.0, B=0.0, C=0.0, D=None, E=None):
        """
        Parameters: *A*, *B*, *C*, *D*, *E*. These must be entered as
        pure numbers without attached dimensions.
        """
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

    def get_yaml(self, out):
        sri = FlowMap([('A', self.A), ('B', self.B), ('C', self.C)])
        if self.D is not None and self.E is not None:
            sri['D'] = self.D
            sri['E'] = self.E
        out['SRI'] = sri


class Lindemann(object):
    """The Lindemann falloff function."""
    def get_yaml(self, out):
        pass


def convert(filename=None, outName=None, text=None):
    import os
    if filename is not None:
        filename = os.path.expanduser(filename)
        base = os.path.basename(filename)
        root, _ = os.path.splitext(base)
        dataset(root)
    if outName is None and _name != 'noname':
        outName = _name + '.yaml'

    open_kw = {'encoding': 'latin-1'} if sys.version_info.major == 3 else {}
    open_mode = 'r' if sys.version_info.major == 3 else 'rU'
    try:
        if filename is not None:
            with open(filename, open_mode, **open_kw) as f:
                code = compile(f.read(), filename, 'exec')
        else:
            code = compile(text, '<string>', 'exec')
        exec(code)
    except SyntaxError as err:
        # Show more context than the default SyntaxError message
        # to help see problems in multi-line statements
        if filename:
            text = open(filename, open_mode).readlines()
        else:
            text = text.split('\n')
        _printerr('%s in "%s" on line %i:\n' % (err.__class__.__name__,
                                                err.filename,
                                                err.lineno))
        _printerr('|  Line |')
        for i in range(max(err.lineno-6, 0),
                       min(err.lineno+3, len(text))):
            _printerr('| % 5i |' % (i+1), text[i].rstrip())
            if i == err.lineno-1:
                _printerr(' '* (err.offset+9) + '^')
        _printerr()
        sys.exit(3)
    except Exception as err:
        import traceback

        if filename:
            text = open(filename, open_mode).readlines()
        else:
            text = text.split('\n')
            filename = '<string>'
        tb = traceback.extract_tb(sys.exc_info()[2])
        lineno = tb[-1][1]
        if tb[-1][0] == filename:
            # Error in input file
            _printerr('%s on line %i of %s:' % (err.__class__.__name__, lineno, filename))
            _printerr(err)
            _printerr('\n| Line |')

            for i in range(max(lineno-6, 0),
                           min(lineno+3, len(text))):
                if i == lineno-1:
                    _printerr('> % 4i >' % (i+1), text[i].rstrip())
                else:
                    _printerr('| % 4i |' % (i+1), text[i].rstrip())
        else:
            # Error in ctml_writer or elsewhere
            traceback.print_exc()

        sys.exit(4)

    # write the YAML file
    emitter = yaml.YAML()
    emitter.width = 70

    for cls in [element, species, thermo, MultiNASA, gas_transport,
                ideal_gas,
                Arrhenius, reaction, three_body_reaction,
                pdep_reaction, falloff_reaction, chemically_activated_reaction,
                chebyshev_reaction, surface_reaction, edge_reaction]:
        emitter.register_class(cls)

    with open(outName, 'w') as dest:
        outputStarted = False

        # todo: header comment
        # todo: motz-wise
        # todo: state
        # todo: fix unit strings to match the new format

        units = FlowMap([])
        if _umass != 'kg':
            units['mass'] = _umass
        if _ulen != 'm':
            units['length'] = _ulen
        if _utime != 's':
            units['time'] = _utime
        if _upres != 'Pa':
            units['pressure'] = _upres
        if _uenergy != 'J':
            units['energy'] = _uenergy
        if _umol != 'kmol':
            units['quantity'] = _umol
        if _ue != 'J/kmol':
            units['activation-energy'] = _ue

        if units:
            unitsMap = BlockMap([('units', units)])
            emitter.dump(unitsMap, dest)
            outputStarted = True

        if _elements:
            elementsMap = BlockMap([('elements', _elements)])
            if outputStarted:
                elementsMap.yaml_set_comment_before_after_key('elements', before='\n')
            outputStarted = True
            emitter.dump(elementsMap, dest)

        if _phases:
            phasesMap = BlockMap([('phases', _phases)])
            if outputStarted:
                phasesMap.yaml_set_comment_before_after_key('phases', before='\n')
            outputStarted = True
            emitter.dump(phasesMap, dest)

        if _species:
            speciesMap = BlockMap([('species', _species)])
            if outputStarted:
                speciesMap.yaml_set_comment_before_after_key('species', before='\n')
            outputStarted = True
            emitter.dump(speciesMap, dest)

        if _reactions:
            reactionsMap = BlockMap([('reactions', _reactions)])
            if outputStarted:
                reactionsMap.yaml_set_comment_before_after_key('reactions', before='\n')
            outputStarted = True
            emitter.dump(reactionsMap, dest)


def main():
    if len(sys.argv) not in (2,3):
        raise ValueError('Incorrect number of command line arguments.')
    convert(*sys.argv[1:])

if __name__ == "__main__":
    main()
