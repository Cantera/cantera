#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains functions for converting Chemkin input files to
Cantera input files (CTI).
"""

import logging
import re
import types

import numpy as np

################################################################################

UNIT_OPTIONS = {'CAL/': 'cal/mol',
                'CAL/MOL': 'cal/mol',
                'CAL/MOLE': 'cal/mol',
                'EVOL': 'eV',
                'EVOLTS': 'eV',
                'JOUL': 'J/mol',
                'JOULES/MOL': 'J/mol',
                'JOULES/MOLE': 'J/mol',
                'KCAL': 'kcal/mol',
                'KCAL/MOL': 'kcal/mol',
                'KCAL/MOLE': 'kcal/mol',
                'KELV': 'K',
                'KELVIN': 'K',
                'KELVINS': 'K',
                'KJOU': 'kJ/mol',
                'KJOULES/MOL': 'kJ/mol',
                'KJOULES/MOLE': 'kJ/mol',
                'MOL': 'mol',
                'MOLE': 'mol',
                'MOLES': 'mol',
                'MOLEC': 'molec',
                'MOLECULES': 'molec'}

ENERGY_UNITS = 'cal/mol'
QUANTITY_UNITS = 'mol'

################################################################################

class ChemkinError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin files. Pass a
    string describing the circumstances that caused the exceptional behavior.
    """
    pass

################################################################################

class Species(object):
    def __init__(self, label):
        self.label = label
        self.thermo = None
        self.transport = None
        self.note = None

    def __str__(self):
        return self.label

    def __repr__(self):
        return 'Species({0!r})'.format(self.label)

    def to_cti(self, indent=0):
        lines = []
        atoms = ' '.join('{0}:{1}'.format(*a) for a in self.composition.iteritems())

        prefix = ' '*(indent+8)

        lines.append('species(name={0!r},'.format(self.label))
        lines.append(prefix + 'atoms={0!r},'.format(atoms))
        if self.thermo:
            lines.append(prefix + 'thermo={0},'.format(self.thermo.to_cti(15+indent)))
        if self.transport:
            lines.append(prefix + 'transport={0},'.format(self.transport.to_cti(14+indent)))
        if self.note:
            lines.append(prefix + 'note={0!r},'.format(self.note))

        lines[-1] = lines[-1][:-1] + ')'
        lines.append('')

        return '\n'.join(lines)

################################################################################

class ThermoModel:
    """
    A base class for thermodynamics models, containing several attributes
    common to all models:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `Tmin`          ``float``           The minimum temperature at which the model is valid, or ``None`` if unknown or undefined
    `Tmax`          ``float``           The maximum temperature at which the model is valid, or ``None`` if unknown or undefined
    `comment`       ``str``             Information about the model (e.g. its source)
    =============== =================== ========================================

    """

    def __init__(self, Tmin=None, Tmax=None, comment=''):
        if Tmin is not None:
            self.Tmin = Tmin
        else:
            self.Tmin = None
        if Tmax is not None:
            self.Tmax = Tmax
        else:
            self.Tmax = None
        self.comment = comment

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ThermoModel object.
        """
        return 'ThermoModel(Tmin={0!r}, Tmax={1!r}, comment="""{2}""")'.format(self.Tmin, self.Tmax, self.comment)

################################################################################

class NASA(ThermoModel):
    """
    A single NASA polynomial for thermodynamic data. The `coeffs` attribute
    stores the seven or nine polynomial coefficients
    :math:`\\mathbf{a} = \\left[a_{-2}\\ a_{-1}\\ a_0\\ a_1\\ a_2\\ a_3\\ a_4\\ a_5\\ a_6 \\right]`
    from which the relevant thermodynamic parameters are evaluated via the
    expressions

    .. math:: \\frac{C_\\mathrm{p}(T)}{R} = a_{-2} T^{-2} + a_{-1} T^{-1} + a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4

    .. math:: \\frac{H(T)}{RT} = - a_{-2} T^{-2} + a_{-1} T^{-1} \\ln T + a_0 + \\frac{1}{2} a_1 T + \\frac{1}{3} a_2 T^2 + \\frac{1}{4} a_3 T^3 + \\frac{1}{5} a_4 T^4 + \\frac{a_5}{T}

    .. math:: \\frac{S(T)}{R} = -\\frac{1}{2} a_{-2} T^{-2} - a_{-1} T^{-1} + a_0 \\ln T + a_1 T + \\frac{1}{2} a_2 T^2 + \\frac{1}{3} a_3 T^3 + \\frac{1}{4} a_4 T^4 + a_6

    The coefficients are stored internally in the nine-coefficient format, even
    when only seven coefficients are provided.
    """

    def __init__(self, coeffs, Tmin=None, Tmax=None, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        coeffs = coeffs or (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        if len(coeffs) == 7:
            self.cm2 = 0.0; self.cm1 = 0.0
            self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
        elif len(coeffs) == 9:
            self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
        else:
            raise ChemkinError('Invalid number of NASA polynomial coefficients; should be 7 or 9.')

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'NASA(Tmin={0!r}, Tmax={1!r}'.format(self.Tmin, self.Tmax)
        if self.cm2 == 0 and self.cm1 == 0:
            string += ', coeffs=[{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g}]'.format(self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
        else:
            string += ', coeffs=[{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g},{7:g},{8:g}]'.format(self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def to_cti(self, indent=0):
        prefix = ' '*indent
        vals = self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6
        vals = ['{0: 15.8E}'.format(i) for i in vals]
        lines = ['NASA([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                 prefix+'     [{0}, {1}, {2},'.format(*vals[0:3]),
                 prefix+'      {0}, {1}, {2},'.format(*vals[3:6]),
                 prefix+'      {0}]),'.format(vals[6])]

        return '\n'.join(lines)

################################################################################

class MultiNASA(ThermoModel):
    """
    A set of thermodynamic parameters given by NASA polynomials. This class
    stores a list of :class:`NASA` objects in the `polynomials`
    attribute. When evaluating a thermodynamic quantity, a polynomial that
    contains the desired temperature within its valid range will be used.
    """

    def __init__(self, polynomials=None, Tmin=0.0, Tmax=0.0, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.polynomials = polynomials or []

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiNASA object.
        """
        string = 'MultiNASA(Tmin={0!r}, Tmax={1!r}'.format(self.Tmin, self.Tmax)
        string += ', polynomials=[{0}]'.format(','.join(['%r' % poly for poly in self.polynomials]))
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def to_cti(self, indent=0):
        prefix = ' '*indent
        lines = []
        for i,p in enumerate(self.polynomials):
            if i == 0:
                lines.append('({0}'.format(p.to_cti(indent+1)))
            elif i != len(self.polynomials)-1:
                lines.append(prefix + ' {0}'.format(p.to_cti(indent+1)))
            else:
                lines.append(prefix + ' {0})'.format(p.to_cti(indent+1)[:-1]))

        return '\n'.join(lines)

################################################################################

class Reaction(object):
    """
    A chemical reaction. The attributes are:

    =================== =========================== ============================
    Attribute           Type                        Description
    =================== =========================== ============================
    `index`             :class:`int`                A unique nonnegative integer index
    `reactants`         :class:`list`               The reactant species (as :class:`Species` objects)
    `products`          :class:`list`               The product species (as :class:`Species` objects)
    `kinetics`          :class:`KineticsModel`      The kinetics model to use for the reaction
    `reversible`        ``bool``                    ``True`` if the reaction is reversible, ``False`` if not
    `transitionState`   :class:`TransitionState`    The transition state
    `thirdBody`         ``bool``                    ``True`` if the reaction if the reaction kinetics imply a third body, ``False`` if not
    `duplicate`         ``bool``                    ``True`` if the reaction is known to be a duplicate, ``False`` if not
    `degeneracy`        :class:`double`             The reaction path degeneracy for the reaction
    `pairs`             ``list``                    Reactant-product pairings to use in converting reaction flux to species flux
    =================== =========================== ============================

    """

    def __init__(self, index=-1, reactants=None, products=None,
                 kinetics=None, reversible=True, transitionState=None,
                 thirdBody=False, duplicate=False, degeneracy=1, pairs=None):
        self.index = index
        self.reactants = reactants
        self.products = products
        self.kinetics = kinetics
        self.reversible = reversible
        self.transitionState = transitionState
        self.thirdBody = thirdBody
        self.duplicate = duplicate
        self.degeneracy = degeneracy
        self.pairs = pairs

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'Reaction('
        if self.index != -1:
            string += 'index={0:d}, '.format(self.index)
        if self.reactants is not None:
            string += 'reactants={0!r}, '.format(self.reactants)
        if self.products is not None:
            string += 'products={0!r}, '.format(self.products)
        if self.kinetics is not None:
            string += 'kinetics={0!r}, '.format(self.kinetics)
        if not self.reversible:
            string += 'reversible={0}, '.format(self.reversible)
        if self.transitionState is not None:
            string += 'transitionState={0!r}, '.format(self.transitionState)
        if self.thirdBody:
            string += 'thirdBody={0}, '.format(self.thirdBody)
        if self.duplicate:
            string += 'duplicate={0}, '.format(self.duplicate)
        if self.degeneracy != 1:
            string += 'degeneracy={0:d}, '.format(self.degeneracy)
        if self.pairs is not None:
            string += 'pairs={0}, '.format(self.pairs)
        string = string[:-2] + ')'
        return string

    def __str__(self):
        """
        Return a string representation of the reaction, in the form 'A + B <=> C + D'.
        """
        arrow = ' <=> '
        if not self.reversible: arrow = ' -> '
        return arrow.join([' + '.join([str(s) for s in self.reactants]), ' + '.join([str(s) for s in self.products])])

    def hasTemplate(self, reactants, products):
        """
        Return ``True`` if the reaction matches the template of `reactants`
        and `products`, which are both lists of :class:`Species` objects, or
        ``False`` if not.
        """
        return ((all([spec in self.reactants for spec in reactants]) and
            all([spec in self.products for spec in products])) or
            (all([spec in self.products for spec in reactants]) and
            all([spec in self.reactants for spec in products])))

    def to_cti(self, indent=0):
        arrow = ' <=> ' if self.reversible else ' => '
        reactantstr = ' + '.join(str(s) for s in self.reactants)
        productstr= ' + '.join(str(s) for s in self.products)

        kinstr = self.kinetics.to_cti(reactantstr, arrow, productstr, indent)


        if self.duplicate:
            k_indent = ' ' * (kinstr.find('(') + 1)
            kinstr = kinstr[:-1] + ",\n{0}options='duplicate')".format(k_indent)

        return kinstr

        if self.thirdBody:
            reactantstr += ' + M'
            productstr += ' + M'
            ctiReactionClass = 'three_body_reaction'
        elif isinstance(self.kinetics, (Lindemann, Troe, Chebyshev)):
            reactantstr += ' (+ M)'
            productstr += ' (+ M)'
            ctiReactionClass = 'falloff_reaction'
        else:
            ctiReactionClass = 'reaction'

        prefix = ' '*(indent+len(self.rxnClass+1))
        reactionstr = reactantstr + arrow + productstr


        if isinstance(self.kinetics, (Arrhenius, ThirdBody)):
            Arates = ' [{0.A}, {0.n}, {0.Ea}]'.format(self.kinetics)
        else:
            Arates = ''

        lines = ['{0}({1!r}{2},'.format(ctiReactionClass, reactionstr, Arates)]

        if self.thirdBody:
            lines.append(prefix + 'efficiencies=')

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

################################################################################

class KineticsModel(object):
    """
    A base class for kinetics models, containing several attributes common to
    all models:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `Tmin`          :class:`Quantity`   The minimum absolute temperature in K at which the model is valid
    `Tmax`          :class:`Quantity`   The maximum absolute temperature in K at which the model is valid
    `Pmin`          :class:`Quantity`   The minimum absolute pressure in Pa at which the model is valid
    `Pmax`          :class:`Quantity`   The maximum absolute pressure in Pa at which the model is valid
    `comment`       :class:`str`        A string containing information about the model (e.g. its source)
    =============== =================== ========================================

    """

    def __init__(self, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        if Tmin is not None:
            self.Tmin = Tmin
        else:
            self.Tmin = None
        if Tmax is not None:
            self.Tmax = Tmax
        else:
            self.Tmax = None
        if Pmin is not None:
            self.Pmin = Pmin
        else:
            self.Pmin = None
        if Pmax is not None:
            self.Pmax = Pmax
        else:
            self.Pmax = None
        self.comment = comment

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KineticsModel object.
        """
        string = self.toPrettyRepr()
        string = re.sub(r'\(\n    ', '(', string)
        string = re.sub(r',\n    ', ', ', string)
        string = re.sub(r',\n\)', ')', string)
        string = re.sub(r' = ', '=', string)
        return string

    def toPrettyRepr(self):
        """
        Return a string representation that can be used to reconstruct the
        KineticsModel object.
        """
        raise NotImplementedError('You must implement this method in your derived class.')

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsModel object.
        """
        return (KineticsModel, (self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isPressureDependent(self):
        """
        Return ``True`` if the kinetics are pressure-dependent or ``False`` if
        they are pressure-independent. This method must be overloaded in the
        derived class.
        """
        raise ChemkinError('Unexpected call to KineticsModel.isPressureDependent(); you should be using a class derived from KineticsModel.')

    def to_cti(self, reactantstr, arrow, productstr):
        raise ChemkinError('to_cti is not implemented for objects of class {0}'.format(self.__class__.__name__))

    def efficiencyString(self):
        if hasattr(self, 'efficiencies'):
            return ' '.join('{0}:{1}'.format(mol, eff)
                            for mol,eff in self.efficiencies.iteritems())
        else:
            return ''

################################################################################

class KineticsData(KineticsModel):
    """
    A kinetics model based around a set of discrete (high-pressure limit)
    rate coefficients at various temperatures. The attributes are:

    =========== =================== ============================================
    Attribute   Type                Description
    =========== =================== ============================================
    `Tdata`     :class:`Quantity`   The temperatures at which the heat capacity data is provided
    `kdata`     :class:`Quantity`   The rate coefficients in SI units at each temperature in `Tdata`
    =========== =================== ============================================

    """

    def __init__(self, Tdata=None, kdata=None, Tmin=None, Tmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.Tdata = Tdata
        self.kdata = kdata

    def toPrettyRepr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = u'KineticsData(\n'
        string += u'    Tdata = {0!r},\n'.format(self.Tdata)
        string += u'    kdata = {0!r},\n'.format(self.kdata)
        if self.Tmin is not None: string += '    Tmin = {0!r},\n'.format(self.Tmin)
        if self.Tmax is not None: string += '    Tmax = {0!r},\n'.format(self.Tmax)
        if self.comment != '': string += '    comment = """{0}""",\n'.format(self.comment)
        return string + u')'

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsData object.
        """
        return (KineticsData, (self.Tdata, self.kdata, self.Tmin, self.Tmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``False`` since KineticsData kinetics are not
        pressure-dependent.
        """
        return False

################################################################################

class Arrhenius(KineticsModel):
    """
    Represent a set of modified Arrhenius kinetics. The kinetic expression has
    the form

    .. math:: k(T) = A \\left( \\frac{T}{T_0} \\right)^n \\exp \\left( - \\frac{E_\\mathrm{a}}{RT} \\right)

    where :math:`A`, :math:`n`, :math:`E_\\mathrm{a}`, and :math:`T_0` are the
    parameters to be set, :math:`T` is absolute temperature, and :math:`R` is
    the gas law constant. The attributes are:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `A`             :class:`Quantity`   The preexponential factor in s^-1, m^3/mol*s, etc.
    `T0`            :class:`Quantity`   The reference temperature in K
    `n`             :class:`Quantity`   The temperature exponent
    `Ea`            :class:`Quantity`   The activation energy in J/mol
    =============== =================== ========================================

    """

    def __init__(self, A=0.0, n=0.0, Ea=0.0, T0=1.0, Tmin=None, Tmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.A = A
        self.T0 = T0
        self.n = n
        self.Ea = Ea

    def toPrettyRepr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = u'Arrhenius(\n'
        string += u'    A = {0!r},\n'.format(self.A)
        string += u'    n = {0!r},\n'.format(self.n)
        string += u'    Ea = {0!r},\n'.format(self.Ea)
        string += u'    T0 = {0!r},\n'.format(self.T0)
        if self.Tmin is not None: string += '    Tmin = {0!r},\n'.format(self.Tmin)
        if self.Tmax is not None: string += '    Tmax = {0!r},\n'.format(self.Tmax)
        if self.comment != '': string += '    comment = """{0}""",\n'.format(self.comment)
        return string + u')'

    def __str__(self):
        """
        Return a string representation that is a bit shorter and prettier than __repr__.
        """
        string = 'Arrhenius(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r})'.format(self.A, self.n, self.Ea, self.T0)
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an Arrhenius object.
        """
        return (Arrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``False`` since Arrhenius kinetics are not pressure-dependent.
        """
        return False

    def rateStr(self):
        return '[{0.A[0]:e}, {0.n}, {0.Ea[0]}]'.format(self)

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstring = reactantstr + arrow + productstr
        return 'reaction({0!r}, {1})'.format(rxnstring, self.rateStr())

################################################################################

class PDepArrhenius(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = A(P) T^{n(P)} \\exp \\left[ \\frac{-E_\\mathrm{a}(P)}{RT} \\right]

    where the modified Arrhenius parameters are stored at a variety of pressures
    and interpolated between on a logarithmic scale. The attributes are:

    =============== ================== ============================================
    Attribute       Type               Description
    =============== ================== ============================================
    `pressures`     :class:`list`      The list of pressures in Pa
    `arrhenius`     :class:`list`      The list of :class:`Arrhenius` objects at each pressure
    `highPlimit`    :class:`Arrhenius` The high (infinite) pressure limiting :class:`Arrhenius` expression
    =============== ================== ============================================

    Note that `highPlimit` is not used in evaluating k(T,P).
    """

    def __init__(self, pressures=None, arrhenius=None, highPlimit=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.pressures = pressures
        self.arrhenius = arrhenius or []
        self.highPlimit = highPlimit or None

    def toPrettyRepr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = u'MultiKinetics(\n'
        string += u'    pressures = {0!r},\n'.format(self.pressures)
        string += u'    arrhenius = [\n'
        for kinetics in self.arrhenius:
            for line in kinetics.toPrettyRepr().splitlines():
                string += u'    {0}\n'.format(line)
        string += u'    ],\n'
        if self.Tmin is not None: string += '    Tmin = {0!r},\n'.format(self.Tmin)
        if self.Tmax is not None: string += '    Tmax = {0!r},\n'.format(self.Tmax)
        if self.Pmin is not None: string += '    Pmin = {0!r},\n'.format(self.Pmin)
        if self.Pmax is not None: string += '    Pmax = {0!r},\n'.format(self.Pmax)
        if self.comment != '': string += '    comment = """{0}""",\n'.format(self.comment)
        return string + u')'

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        PDepArrhenius object.
        """
        string = 'PDepArrhenius(\n pressures={0!r},\n arrhenius=[\n  {1}]'.format(self.pressures, ',\n  '.join([repr(arrh) for arrh in self.arrhenius]))
        if self.highPlimit is not None: string += ",\n highPlimit={0!r}".format(self.highPlimit)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ',\n comment="""{0}"""'.format(self.comment)
        string += '\n)'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a PDepArrhenius object.
        """
        return (PDepArrhenius, (self.pressures, self.arrhenius, self.highPlimit, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``True`` since PDepArrhenius kinetics are pressure-dependent.
        """
        return True

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstring = reactantstr + arrow + productstr
        lines = ['pdep_arrhenius({0!r},'.format(rxnstring)]
        prefix = ' '*(indent+15)
        template = '[({0}, {1!r}), {2.A[0]:e}, {2.n}, {2.Ea[0]}],'
        for pressure,arrhenius in zip(self.pressures[0], self.arrhenius):
            lines.append(prefix + template.format(pressure,
                                                  self.pressures[1],
                                                  arrhenius))
        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

class Chebyshev(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: \\log k(T,P) = \\sum_{t=1}^{N_T} \\sum_{p=1}^{N_P} \\alpha_{tp} \\phi_t(\\tilde{T}) \\phi_p(\\tilde{P})

    where :math:`\\alpha_{tp}` is a constant, :math:`\\phi_n(x)` is the
    Chebyshev polynomial of degree :math:`n` evaluated at :math:`x`, and

    .. math:: \\tilde{T} \\equiv \\frac{2T^{-1} - T_\\mathrm{min}^{-1} - T_\\mathrm{max}^{-1}}{T_\\mathrm{max}^{-1} - T_\\mathrm{min}^{-1}}

    .. math:: \\tilde{P} \\equiv \\frac{2 \\log P - \\log P_\\mathrm{min} - \\log P_\\mathrm{max}}{\\log P_\\mathrm{max} - \\log P_\\mathrm{min}}

    are reduced temperature and reduced pressures designed to map the ranges
    :math:`(T_\\mathrm{min}, T_\\mathrm{max})` and
    :math:`(P_\\mathrm{min}, P_\\mathrm{max})` to :math:`(-1, 1)`.
    The attributes are:

    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `coeffs`        :class:`list`   Matrix of Chebyshev coefficients
    `kunits`        ``str``         The units of the generated k(T, P) values
    `degreeT`       :class:`int`    The number of terms in the inverse temperature direction
    `degreeP`       :class:`int`    The number of terms in the log pressure direction
    =============== =============== ============================================

    """

    def __init__(self, coeffs=None, kunits='', Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        if coeffs is not None:
            self.coeffs = np.array(coeffs, np.float64)
            self.degreeT = self.coeffs.shape[0]
            self.degreeP = self.coeffs.shape[1]
        else:
            self.coeffs = None
            self.degreeT = 0
            self.degreeP = 0
        self.kunits = kunits

    def toPrettyRepr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = u'Chebyshev(\n'
        string += u'    coeffs = [\n'
        for i in range(self.degreeT):
            string += u'        [{0}]'.format(','.join(['{0:g}'.format(self.coeffs[i,j]) for j in range(self.degreeP)]))
        string += u'    ],\n'
        if self.Tmin is not None: string += '    Tmin = {0!r},\n'.format(self.Tmin)
        if self.Tmax is not None: string += '    Tmax = {0!r},\n'.format(self.Tmax)
        if self.Pmin is not None: string += '    Pmin = {0!r},\n'.format(self.Pmin)
        if self.Pmax is not None: string += '    Pmax = {0!r},\n'.format(self.Pmax)
        if self.comment != '': string += '    comment = """{0}""",\n'.format(self.comment)
        return string + u')'

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Chebyshev object.
        """
        coeffs = '['
        for i in range(self.degreeT):
            if i > 0: coeffs += ', '
            coeffs += '[{0}]'.format(','.join(['{0:g}'.format(self.coeffs[i,j]) for j in range(self.degreeP)]))
        coeffs += ']'

        string = 'Chebyshev(coeffs={0}'.format(coeffs)
        if self.kunits != '': string += ', kunits="{0}"'.format(self.kunits)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Chebyshev object.
        """
        return (Chebyshev, (self.coeffs, self.kunits, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``True`` since Chebyshev polynomial kinetics are
        pressure-dependent.
        """
        return True

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + ' (+ M)' + arrow + productstr + ' (+ M)'
        prefix = ' '*(indent+19)
        lines = ['chebyshev_reaction({0!r},'.format(rxnstr),
                 prefix + 'Tmin={0.Tmin}, Tmax={0.Tmax},'.format(self),
                 prefix + 'Pmin={0.Pmin}, Pmax={0.Pmax},'.format(self)]
        for i in range(self.degreeT):
            coeffline = ', '.join('{0: 12.5e}'.format(self.coeffs[i,j]) for j in range(self.degreeP))
            if i == 0:
                lines.append(prefix + 'coeffs=[[{0}],'.format(coeffline))
            else:
                lines.append(prefix + '        [{0}],'.format(coeffline))

        lines[-1] = lines[-1][:-1] + '])'
        return '\n'.join(lines)

################################################################################

class ThirdBody(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k(T) [\\ce{M}]

    where :math:`k(T)` is an Arrhenius expression and
    :math:`[\\ce{M}] \\approx P/RT` is the concentration of the third body
    (i.e. the bath gas). A collision efficiency can be used to further correct
    the value of :math:`k(T,P)`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusHigh=None, efficiencies=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrheniusHigh = arrheniusHigh
        self.efficiencies = {}
        if efficiencies is not None:
            for mol, eff in efficiencies.iteritems():
                self.efficiencies[mol] = eff

    def toPrettyRepr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = u'ThirdBody(\n'

        lines = self.arrheniusHigh.toPrettyRepr().splitlines()
        string += u'    arrheniusHigh = {0}\n'.format(lines[0])
        for line in lines[1:-1]:
            string += u'    {0}\n'.format(line)
        string += u'    ),\n'

        if len(self.efficiencies) > 0:
            string += u'    efficiencies = {\n'
            for species in sorted(self.efficiencies):
                string += u'        "{0}": {1:g},\n'.format(species, self.efficiencies[species])
            string += u'    },\n'

        if self.Tmin is not None: string += '    Tmin = {0!r},\n'.format(self.Tmin)
        if self.Tmax is not None: string += '    Tmax = {0!r},\n'.format(self.Tmax)
        if self.Pmin is not None: string += '    Pmin = {0!r},\n'.format(self.Pmin)
        if self.Pmax is not None: string += '    Pmax = {0!r},\n'.format(self.Pmax)
        if self.comment != '': string += '    comment = """{0}""",\n'.format(self.comment)
        return string + u')'

    def __reduce__(self):
        """
        A helper function used when pickling a ThirdBody object.
        """
        return (ThirdBody, (self.arrheniusHigh, self.efficiencies, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``True`` since third-body kinetics are pressure-dependent.
        """
        return True

    def getColliderEfficiency(self, collider):
        """
        Return the collider efficiency for the specified `collider`, which can
        take one of two forms:

        * A single collider species. If the collider exists in the in the set
          of efficiencies, its efficiency will be returned. If not, an
          efficiency of unity will be returned.

        * A ``dict`` mapping collider species to mole fractions. The overall
          efficiency will be a weighted sum of the efficiencies of the collider
          species, using the mole fractions as the weights. Collider species not
          present in the set of efficiencies will be assumed to have an
          efficiency of unity.

        If collider is ``None`` or otherwise invalid, an efficiency of unity
        will be returned.
        """
        if isinstance(collider, dict):
            # Assume collider is a dict mapping species to weights
            efficiency = 0.0
            for spec, frac in collider.iteritems():
                try:
                    eff = self.efficiencies[spec]
                except KeyError:
                    eff = 1.0
                efficiency += eff * frac
            efficiency /= sum(collider.values())
        else:
            # Assume collider is a single species
            try:
                efficiency = self.efficiencies[collider]
            except KeyError:
                efficiency = 1.0

        return efficiency

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + ' + M' + arrow + productstr + ' + M'
        prefix = ' '*(indent + 20)
        lines = ['three_body_reaction({0!r}, {1},'.format(rxnstr, self.arrheniusHigh.rateStr())]
        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

class Lindemann(ThirdBody):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k_\\infty(T) \\left[ \\frac{P_\\mathrm{r}}{1 + P_\\mathrm{r}} \\right] F

    where

    .. math::

        P_\\mathrm{r} &= \\frac{k_0(T)}{k_\\infty(T)} [\\ce{M}]

        k_0(T) &= A_0 T^{n_0} \\exp \\left( - \\frac{E_0}{RT} \\right)

        k_\\infty(T) &= A_\\infty T^{n_\\infty} \\exp \\left( - \\frac{E_\\infty}{RT} \\right)

    and :math:`[\\ce{M}] \\approx P/RT` is the concentration of the
    bath gas. The Arrhenius expressions :math:`k_0(T)` and :math:`k_\\infty(T)`
    represent the low-pressure and high-pressure limit kinetics, respectively.
    The former is necessarily one reaction order higher than the latter. For
    the Lindemann model, :math:`F = 1`. A collision efficiency can be used to
    further correct the value of :math:`k(T,P)`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusLow=None, arrheniusHigh=None, efficiencies=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        ThirdBody.__init__(self, arrheniusHigh=arrheniusHigh, efficiencies=efficiencies, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrheniusLow = arrheniusLow

    def toPrettyRepr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = u'Lindemann(\n'

        lines = self.arrheniusHigh.toPrettyRepr().splitlines()
        string += u'    arrheniusHigh = {0}\n'.format(lines[0])
        for line in lines[1:-1]:
            string += u'    {0}\n'.format(line)
        string += u'    ),\n'

        lines = self.arrheniusLow.toPrettyRepr().splitlines()
        string += u'    arrheniusLow = {0}\n'.format(lines[0])
        for line in lines[1:-1]:
            string += u'    {0}\n'.format(line)
        string += u'    ),\n'

        if len(self.efficiencies) > 0:
            string += u'    efficiencies = {\n'
            for species in sorted(self.efficiencies):
                string += u'        "{0}": {1:g},\n'.format(species, self.efficiencies[species])
            string += u'    },\n'

        if self.Tmin is not None: string += '    Tmin = {0!r},\n'.format(self.Tmin)
        if self.Tmax is not None: string += '    Tmax = {0!r},\n'.format(self.Tmax)
        if self.Pmin is not None: string += '    Pmin = {0!r},\n'.format(self.Pmin)
        if self.Pmax is not None: string += '    Pmax = {0!r},\n'.format(self.Pmax)
        if self.comment != '': string += '    comment = """{0}""",\n'.format(self.comment)
        return string + u')'

    def __reduce__(self):
        """
        A helper function used when pickling a Lindemann object.
        """
        return (Lindemann, (self.arrheniusLow, self.arrheniusHigh, self.efficiencies, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + ' (+ M)' + arrow + productstr + ' (+ M)'
        prefix = ' '*(indent + 17)
        lines = ['falloff_reaction({0!r},'.format(rxnstr)]
        lines.append(prefix + 'kf={0},'.format(self.arrheniusHigh.rateStr()))
        lines.append(prefix + 'kf0={0},'.format(self.arrheniusLow.rateStr()))
        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

class Troe(Lindemann):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k_\\infty(T) \\left[ \\frac{P_\\mathrm{r}}{1 + P_\\mathrm{r}} \\right] F

    where

    .. math::

        P_\\mathrm{r} &= \\frac{k_0(T)}{k_\\infty(T)} [\\ce{M}]

        k_0(T) &= A_0 T^{n_0} \\exp \\left( - \\frac{E_0}{RT} \\right)

        k_\\infty(T) &= A_\\infty T^{n_\\infty} \\exp \\left( - \\frac{E_\\infty}{RT} \\right)

    and :math:`[\\ce{M}] \\approx P/RT` is the concentration of the
    bath gas. The Arrhenius expressions :math:`k_0(T)` and :math:`k_\\infty(T)`
    represent the low-pressure and high-pressure limit kinetics, respectively.
    The former is necessarily one reaction order higher than the latter. A
    collision efficiency can be used to further correct the value of
    :math:`k(T,P)`.

    For the Troe model the parameter :math:`F` is computed via

    .. math::

        \\log F &= \\left\\{1 + \\left[ \\frac{\\log P_\\mathrm{r} + c}{n - d (\\log P_\\mathrm{r} + c)} \\right]^2 \\right\\}^{-1} \\log F_\\mathrm{cent}

        c &= -0.4 - 0.67 \\log F_\\mathrm{cent}

        n &= 0.75 - 1.27 \\log F_\\mathrm{cent}

        d &= 0.14

        F_\\mathrm{cent} &= (1 - \\alpha) \\exp \\left( -T/T_3 \\right) + \\alpha \\exp \\left( -T/T_1 \\right) + \\exp \\left( -T_2/T \\right)

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    `alpha`         :class:`Quantity`       The :math:`\\alpha` parameter
    `T1`            :class:`Quantity`       The :math:`T_1` parameter
    `T2`            :class:`Quantity`       The :math:`T_2` parameter
    `T3`            :class:`Quantity`       The :math:`T_3` parameter
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusLow=None, arrheniusHigh=None, efficiencies=None, alpha=0.0, T3=0.0, T1=0.0, T2=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        Lindemann.__init__(self, arrheniusLow=arrheniusLow, arrheniusHigh=arrheniusHigh, efficiencies=efficiencies, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.alpha = alpha
        self.T3 = T3
        self.T1 = T1
        if T2 is not None:
            self.T2 = T2
        else:
            self.T2 = None

    def toPrettyRepr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = u'Troe(\n'

        lines = self.arrheniusHigh.toPrettyRepr().splitlines()
        string += u'    arrheniusHigh = {0}\n'.format(lines[0])
        for line in lines[1:-1]:
            string += u'    {0}\n'.format(line)
        string += u'    ),\n'

        lines = self.arrheniusLow.toPrettyRepr().splitlines()
        string += u'    arrheniusLow = {0}\n'.format(lines[0])
        for line in lines[1:-1]:
            string += u'    {0}\n'.format(line)
        string += u'    ),\n'

        string += u'    alpha = {0!r},\n'.format(self.alpha)
        string += u'    T3 = {0!r},\n'.format(self.T3)
        string += u'    T1 = {0!r},\n'.format(self.T1)
        if self.T2 is not None: string += u'    T2 = {0!r},\n'.format(self.T2)

        if len(self.efficiencies) > 0:
            string += u'    efficiencies = {\n'
            for molecule in sorted(self.efficiencies):
                string += u'        "{0}": {1:g},\n'.format(molecule, self.efficiencies[molecule])
            string += u'    },\n'

        if self.Tmin is not None: string += '    Tmin = {0!r},\n'.format(self.Tmin)
        if self.Tmax is not None: string += '    Tmax = {0!r},\n'.format(self.Tmax)
        if self.Pmin is not None: string += '    Pmin = {0!r},\n'.format(self.Pmin)
        if self.Pmax is not None: string += '    Pmax = {0!r},\n'.format(self.Pmax)
        if self.comment != '': string += '    comment = """{0}""",\n'.format(self.comment)
        return string + u')'

    def __reduce__(self):
        """
        A helper function used when pickling a Troe object.
        """
        return (Troe, (self.arrheniusLow, self.arrheniusHigh, self.efficiencies, self.alpha, self.T3, self.T1, self.T2, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + ' (+ M)' + arrow + productstr + ' (+ M)'
        prefix = ' '*17
        lines = ['falloff_reaction({0!r},'.format(rxnstr),
                 prefix + 'kf={0},'.format(self.arrheniusHigh.rateStr()),
                 prefix + 'kf0={0},'.format(self.arrheniusLow.rateStr())]

        if self.T2:
            troeArgs = 'A={0.alpha[0]}, T3={0.T3[0]}, T1={0.T1[0]}, T2={0.T2[0]}'.format(self)
        else:
            troeArgs = 'A={0.alpha[0]}, T3={0.T3[0]}, T1={0.T1[0]}'.format(self)
        lines.append(prefix + 'falloff=Troe({0}),'.format(troeArgs))

        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))

        # replace trailing comma
        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

class TransportData(object):
    geometryFlags = ['atom', 'linear', 'nonlinear']

    def __init__(self, label, geometry, wellDepth, collisionDiameter,
                 dipoleMoment, polarizability, zRot, comment=None):

        assert isinstance(label, types.StringTypes)
        assert int(geometry) in (0,1,2)

        self.label = label
        self.geometry = self.geometryFlags[int(geometry)]
        self.wellDepth = float(wellDepth)
        self.collisionDiameter = float(collisionDiameter)
        self.dipoleMoment = float(dipoleMoment)
        self.polarizability = float(polarizability)
        self.zRot = float(zRot)
        self.comment = comment or '' # @todo: include this in the CTI

    def __repr__(self):
        return ('TransportData({label!r}, {geometry!r}, {wellDepth!r}, '
                '{collisionDiameter!r}, {dipoleMoment!r}, {polarizability!r}, '
                '{zRot!r}, {comment!r})').format(**self.__dict__)

    def to_cti(self, indent=0):
        prefix = ' '*(indent+18)
        lines = ['gas_transport(geom={0!r},'.format(self.geometry),
                 prefix+'diam={0},'.format(self.collisionDiameter),
                 prefix+'well_depth={0},'.format(self.wellDepth)]
        if self.dipoleMoment:
            lines.append(prefix+'dipole={0},'.format(self.dipoleMoment))
        if self.polarizability:
            lines.append(prefix+'polar={0},'.format(self.polarizability))
        if self.zRot:
            lines.append(prefix+'rot_relax={0},'.format(self.zRot))

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

def readThermoEntry(entry):
    """
    Read a thermodynamics `entry` for one species in a Chemkin file. Returns
    the label of the species, the thermodynamics model as a :class:`MultiNASA`
    object and the elemental composition of the species.
    """
    lines = entry.splitlines()
    identifier = lines[0][0:24].split()
    species = identifier[0].strip()

    if len(identifier) > 1:
        note = ''.join(identifier[1:]).strip()
    else:
        note = ''

    # Extract the NASA polynomial coefficients
    # Remember that the high-T polynomial comes first!
    try:
        Tmin = float(lines[0][45:55].strip())
        Tmax = float(lines[0][55:65].strip())
        Tint = float(lines[0][65:75].strip())

        a0_high = float(lines[1][0:15].strip())
        a1_high = float(lines[1][15:30].strip())
        a2_high = float(lines[1][30:45].strip())
        a3_high = float(lines[1][45:60].strip())
        a4_high = float(lines[1][60:75].strip())

        a5_high = float(lines[2][0:15].strip())
        a6_high = float(lines[2][15:30].strip())
        a0_low = float(lines[2][30:45].strip())
        a1_low = float(lines[2][45:60].strip())
        a2_low = float(lines[2][60:75].strip())

        a3_low = float(lines[3][0:15].strip())
        a4_low = float(lines[3][15:30].strip())
        a5_low = float(lines[3][30:45].strip())
        a6_low = float(lines[3][45:60].strip())
    except (IndexError, ValueError):
        raise ChemkinError('Error while reading thermo entry for species {0}'.format(species))

    elements = lines[0][24:44]
    composition = {}
    for i in range(4):
        symbol = elements[5*i:5*i+2].strip()
        count = elements[5*i+2:5*i+5].strip()
        if not symbol:
            continue
        try:
            count = int(float(count))
            if count:
                composition[symbol.capitalize()] = count
        except ValueError:
            pass

    # Construct and return the thermodynamics model
    thermo = MultiNASA(
        polynomials = [
            NASA(Tmin=(Tmin,"K"), Tmax=(Tint,"K"), coeffs=[a0_low, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low]),
            NASA(Tmin=(Tint,"K"), Tmax=(Tmax,"K"), coeffs=[a0_high, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high])
        ],
        Tmin = (Tmin,"K"),
        Tmax = (Tmax,"K"),
    )

    return species, thermo, composition, note

################################################################################

def readKineticsEntry(entry, speciesDict, energyUnits, moleculeUnits):
    """
    Read a kinetics `entry` for a single reaction as loaded from a Chemkin
    file. The associated mapping of labels to species `speciesDict` should also
    be provided. Returns a :class:`Reaction` object with the reaction and its
    associated kinetics.
    """

    lines = entry.strip().splitlines()

    # The first line contains the reaction equation and a set of modified Arrhenius parameters
    tokens = lines[0].split()
    A = float(tokens[-3])
    n = float(tokens[-2])
    Ea = float(tokens[-1])
    reaction = ''.join(tokens[:-3])
    thirdBody = False

    # Split the reaction equation into reactants and products
    if '<=>' in reaction:
        reversible = True
        reactants, products = reaction.split('<=>')
    elif '=>' in reaction:
        reversible = False
        reactants, products = reaction.split('=>')
    elif '=' in reaction:
        reversible = True
        reactants, products = reaction.split('=')
    else:
        raise ChemkinError("Failed to find reactant/product delimiter in reaction string.")

    if '(+M)' in reactants: reactants = reactants.replace('(+M)','')
    if '(+m)' in reactants: reactants = reactants.replace('(+m)','')
    if '(+M)' in products:  products = products.replace('(+M)','')
    if '(+m)' in products:  products = products.replace('(+m)','')

    # Create a new Reaction object for this reaction
    reaction = Reaction(reactants=[], products=[], reversible=reversible)

    # Convert the reactants and products to Species objects using the speciesDict
    for reactant in reactants.split('+'):
        reactant = reactant.strip()
        stoichiometry = 1
        if reactant[0].isdigit():
            # This allows for reactions to be of the form 2A=B+C instead of A+A=B+C
            # The implementation below assumes an integer between 0 and 9, inclusive
            stoichiometry = int(reactant[0])
            reactant = reactant[1:]
        if reactant == 'M' or reactant == 'm':
            thirdBody = True
        elif reactant not in speciesDict:
            raise ChemkinError('Unexpected reactant "{0}" in reaction {1}.'.format(reactant, reaction))
        else:
            for i in range(stoichiometry):
                reaction.reactants.append(speciesDict[reactant])
    for product in products.split('+'):
        product = product.strip()
        stoichiometry = 1
        if product[0].isdigit():
            # This allows for reactions to be of the form A+B=2C instead of A+B=C+C
            # The implementation below assumes an integer between 0 and 9, inclusive
            stoichiometry = int(product[0])
            product = product[1:]
        if product.upper() == 'M' or product == 'm':
            pass
        elif product not in speciesDict:
            raise ChemkinError('Unexpected product "{0}" in reaction {1}.'.format(product, reaction))
        else:
            for i in range(stoichiometry):
                reaction.products.append(speciesDict[product])

    # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
    # This assumes elementary kinetics for all reactions
    if len(reaction.reactants) + (1 if thirdBody else 0) == 3:
        kunits = "cm^6/(mol^2*s)"
        klow_units = "cm^9/(mol^3*s)"
    elif len(reaction.reactants) + (1 if thirdBody else 0) == 2:
        kunits = "cm^3/(mol*s)"
        klow_units = "cm^6/(mol^2*s)"
    elif len(reaction.reactants) + (1 if thirdBody else 0) == 1:
        kunits = "s^-1"
        klow_units = "cm^3/(mol*s)"
    else:
        raise ChemkinError('Invalid number of reactant species for reaction {0}.'.format(reaction))

    # The rest of the first line contains the high-P limit Arrhenius parameters (if available)
    #tokens = lines[0][52:].split()
    tokens = lines[0].split()[1:]
    arrheniusHigh = Arrhenius(
        A = (A,kunits),
        n = n,
        Ea = (Ea, energyUnits),
        T0 = (1,"K"),
    )

    if len(lines) == 1:
        # If there's only one line then we know to use the high-P limit kinetics as-is
        reaction.kinetics = arrheniusHigh
    else:
        # There's more kinetics information to be read
        arrheniusLow = None
        troe = None
        chebyshev = None
        pdepArrhenius = None
        efficiencies = {}
        chebyshevCoeffs = []

        # Note that the subsequent lines could be in any order
        for line in lines[1:]:

            tokens = line.split('/')
            if 'DUP' in line or 'dup' in line:
                # Duplicate reaction
                reaction.duplicate = True

            elif 'LOW' in line or 'low' in line:
                # Low-pressure-limit Arrhenius parameters
                tokens = tokens[1].split()
                arrheniusLow = Arrhenius(
                    A = (float(tokens[0].strip()),klow_units),
                    n = float(tokens[1].strip()),
                    Ea = (float(tokens[2].strip()),"kcal/mol"),
                    T0 = (1,"K"),
                )

            elif 'TROE' in line or 'troe' in line:
                # Troe falloff parameters
                tokens = tokens[1].split()
                alpha = float(tokens[0].strip())
                T3 = float(tokens[1].strip())
                T1 = float(tokens[2].strip())
                try:
                    T2 = float(tokens[3].strip())
                except (IndexError, ValueError):
                    T2 = None

                troe = Troe(
                    alpha = (alpha,''),
                    T3 = (T3,"K"),
                    T1 = (T1,"K"),
                    T2 = (T2,"K") if T2 is not None else None,
                )

            elif 'CHEB' in line or 'cheb' in line:
                # Chebyshev parameters
                if chebyshev is None:
                    chebyshev = Chebyshev()
                tokens = [t.strip() for t in tokens]
                if 'TCHEB' in line:
                    index = tokens.index('TCHEB')
                    tokens2 = tokens[index+1].split()
                    chebyshev.Tmin = float(tokens2[0].strip())
                    chebyshev.Tmax = float(tokens2[1].strip())
                if 'PCHEB' in line:
                    index = tokens.index('PCHEB')
                    tokens2 = tokens[index+1].split()
                    chebyshev.Pmin = (float(tokens2[0].strip()), 'atm')
                    chebyshev.Pmax = (float(tokens2[1].strip()), 'atm')
                if 'TCHEB' in line or 'PCHEB' in line:
                    pass
                elif chebyshev.degreeT == 0 or chebyshev.degreeP == 0:
                    tokens2 = tokens[1].split()
                    chebyshev.degreeT = int(float(tokens2[0].strip()))
                    chebyshev.degreeP = int(float(tokens2[1].strip()))
                    chebyshev.coeffs = np.zeros((chebyshev.degreeT,chebyshev.degreeP), np.float64)
                else:
                    tokens2 = tokens[1].split()
                    chebyshevCoeffs.extend([float(t.strip()) for t in tokens2])

            elif 'PLOG' in line or 'plog' in line:
                # Pressure-dependent Arrhenius parameters
                if pdepArrhenius is None:
                    pdepArrhenius = []
                tokens = tokens[1].split()
                pdepArrhenius.append([float(tokens[0].strip()), Arrhenius(
                    A = (float(tokens[1].strip()),kunits),
                    n = float(tokens[2].strip()),
                    Ea = (float(tokens[3].strip()),"kcal/mol"),
                    T0 = (1,"K"),
                )])

            else:
                # Assume a list of collider efficiencies
                for collider, efficiency in zip(tokens[0::2], tokens[1::2]):
                    efficiencies[collider.strip()] = float(efficiency.strip())

        # Decide which kinetics to keep and store them on the reaction object
        # Only one of these should be true at a time!
        if chebyshev is not None:
            if chebyshev.Tmin is None or chebyshev.Tmax is None:
                raise ChemkinError('Missing TCHEB line for reaction {0}'.format(reaction))
            if chebyshev.Pmin is None or chebyshev.Pmax is None:
                raise ChemkinError('Missing PCHEB line for reaction {0}'.format(reaction))
            index = 0
            for t in range(chebyshev.degreeT):
                for p in range(chebyshev.degreeP):
                    chebyshev.coeffs[t,p] = chebyshevCoeffs[index]
                    index += 1
            reaction.kinetics = chebyshev
        elif pdepArrhenius is not None:
            reaction.kinetics = PDepArrhenius(
                pressures = ([P for P, arrh in pdepArrhenius],"atm"),
                arrhenius = [arrh for P, arrh in pdepArrhenius],
            )
        elif troe is not None:
            troe.arrheniusHigh = arrheniusHigh
            troe.arrheniusLow = arrheniusLow
            troe.efficiencies = efficiencies
            reaction.kinetics = troe
        elif arrheniusLow is not None:
            reaction.kinetics = Lindemann(arrheniusHigh=arrheniusHigh, arrheniusLow=arrheniusLow)
            reaction.kinetics.efficiencies = efficiencies
        elif thirdBody:
            reaction.kinetics = ThirdBody(arrheniusHigh=arrheniusHigh)
            reaction.kinetics.efficiencies = efficiencies
        elif reaction.duplicate:
            reaction.kinetics = arrheniusHigh
        else:
            raise ChemkinError('Unable to determine pressure-dependent kinetics for reaction {0}.'.format(reaction))

    return reaction

################################################################################

def loadChemkinFile(path):
    """
    Load a Chemkin input file to `path` on disk, returning lists of the species
    and reactions in the Chemkin file.
    """

    speciesList = []; speciesDict = {}
    reactionList = []
    transportLines = []

    def removeCommentFromLine(line):
        if '!' in line:
            index = line.index('!')
            comment = line[index+1:-1]
            line = line[0:index] + '\n'
            return line, comment
        else:
            comment = ''
            return line, comment

    with open(path, 'r') as f:
        line = f.readline()
        while line != '':
            line = removeCommentFromLine(line)[0]
            line = line.strip()
            tokens = line.split()

            if 'SPECIES' in line:
                # List of species identifiers
                index = tokens.index('SPECIES')
                tokens = tokens[index+1:]
                while 'END' not in tokens:
                    line = f.readline()
                    line = removeCommentFromLine(line)[0]
                    line = line.strip()
                    tokens.extend(line.split())

                for token in tokens:
                    if token == 'END':
                        break
                    if token in speciesDict:
                        species = speciesDict[token]
                    else:
                        species = Species(label=token)
                        speciesDict[token] = species
                    speciesList.append(species)

            elif 'THERM' in line:
                # List of thermodynamics (hopefully one per species!)
                line = f.readline()
                thermo = ''
                while line != '' and 'END' not in line:
                    line = removeCommentFromLine(line)[0]
                    if len(line) >= 80:
                        if line[79] in ['1', '2', '3', '4']:
                            thermo += line
                            if line[79] == '4':
                                label, thermo, comp, note = readThermoEntry(thermo)
                                try:
                                    speciesDict[label].thermo = thermo
                                    speciesDict[label].composition = comp
                                    speciesDict[label].note = note
                                except KeyError:
                                    logging.warning('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))
                                thermo = ''
                    line = f.readline()

            elif 'REACTIONS' in line:
                # Reactions section
                energyUnits = 'CAL/MOL'
                moleculeUnits = 'MOLES'
                try:
                    energyUnits = tokens[1]
                    moleculeUnits = tokens[2]
                except IndexError:
                    pass

                ENERGY_UNITS = UNIT_OPTIONS[energyUnits]
                QUANTITY_UNITS = UNIT_OPTIONS[moleculeUnits]

                kineticsList = []
                commentsList = []
                kinetics = ''
                comments = ''

                line = f.readline()
                while line != '' and 'END' not in line:

                    lineStartsWithComment = line.startswith('!')
                    line, comment = removeCommentFromLine(line)
                    line = line.strip(); comment = comment.strip()

                    if 'rev' in line or 'REV' in line:
                        # can no longer name reactants rev...
                        line = f.readline()

                    if '=' in line and not lineStartsWithComment:
                        # Finish previous record
                        kineticsList.append(kinetics)
                        commentsList.append(comments)
                        kinetics = ''
                        comments = ''

                    if line: kinetics += line + '\n'
                    if comment: comments += comment + '\n'

                    line = f.readline()

                # Don't forget the last reaction!
                if kinetics.strip() != '':
                    kineticsList.append(kinetics)
                    commentsList.append(comments)

                if kineticsList[0] == '' and commentsList[-1] == '':
                    # True for Chemkin files generated from RMG-Py
                    kineticsList.pop(0)
                    commentsList.pop(-1)
                elif kineticsList[0] == '' and commentsList[0] == '':
                    # True for Chemkin files generated from RMG-Java
                    kineticsList.pop(0)
                    commentsList.pop(0)
                else:
                    # In reality, comments can occur anywhere in the Chemkin
                    # file (e.g. either or both of before and after the
                    # reaction equation)
                    # If we can't tell what semantics we are using, then just
                    # throw the comments away
                    # (This is better than failing to load the Chemkin file at
                    # all, which would likely occur otherwise)
                    if kineticsList[0] == '':
                        kineticsList.pop(0)
                    if len(kineticsList) != len(commentsList):
                        commentsList = ['' for kinetics in kineticsList]

                for kinetics, comments in zip(kineticsList, commentsList):
                    reaction = readKineticsEntry(kinetics, speciesDict, energyUnits, moleculeUnits)
                    reactionList.append(reaction)

            elif 'TRAN' in line:
                line = f.readline()
                while 'END' not in line:
                    transportLines.append(line)

            line = f.readline()

    # Check for marked (and unmarked!) duplicate reactions
    # Raise exception for unmarked duplicate reactions
    for index1 in range(len(reactionList)):
        reaction1 = reactionList[index1]
        for index2 in range(index1+1, len(reactionList)):
            reaction2 = reactionList[index2]
            if reaction1.reactants == reaction2.reactants and reaction1.products == reaction2.products:
                if reaction1.duplicate and reaction2.duplicate:
                    pass
                elif reaction1.kinetics.isPressureDependent() == reaction2.kinetics.isPressureDependent():
                    # If both reactions are pressure-independent or both are pressure-dependent, then they need duplicate tags
                    # Chemkin treates pdep and non-pdep reactions as different, so those are okay
                    raise ChemkinError('Encountered unmarked duplicate reaction {0}.'.format(reaction1))

    index = 0
    for reaction in reactionList:
        index += 1
        reaction.index = index

    if transportLines:
        parseTransportData(transportLines, speciesList)

    return speciesList, reactionList

################################################################################

def parseTransportData(lines, speciesList):
    """
    Parse the Chemkin-format transport data in ``lines`` (a list of strings)
    and add that transport data to the species in ``speciesList``.
    """
    speciesDict = dict((species.label, species) for species in speciesList)
    for line in lines:
        line = line.strip()
        if not line or line.startswith('!'):
            continue

        data = line.split()
        if len(data) < 7:
            raise ChemkinError('Unable to parse transport data: not enough parameters')
        if len(data) >= 8:
            # comment may contain spaces. Rejoin into a single field.
            comment = ''.join(data[7:]).lstrip('!')
            data = data[:7] + [comment]

        speciesName = data[0]
        if speciesName in speciesDict:
            speciesDict[speciesName].transport = TransportData(*data)

################################################################################

def writeCTI(species,
             reactions=None,
             header=None,
             name='gas',
             transportModel='Mix',
             outName='mech.cti'):

    delimiterLine = '#' + '-'*79
    haveTransport = True
    speciesNameLength = 1
    elements = set()
    for s in species:
        if not s.transport:
            haveTransport = False
        elements.update(s.composition)
        speciesNameLength = max(speciesNameLength, len(s.label))

    speciesNames = ['']
    for i,s in enumerate(species):
        if i and not i % 5:
            speciesNames.append(' '*21)
        speciesNames[-1] += '{0:{1}s}'.format(s.label, speciesNameLength+2)

    speciesNames = '\n'.join(speciesNames).strip()

    lines = []
    if header:
        lines.extend(header)

    # Write the gas definition
    lines.append("units(length='cm', time='s', quantity={0!r}, act_energy={1!r})".format(QUANTITY_UNITS, ENERGY_UNITS))
    lines.append('')
    lines.append('ideal_gas(name={0!r},'.format(name))
    lines.append('          elements="{0}",'.format(' '.join(elements)))
    lines.append('          species="""{0}""",'.format(speciesNames))
    if reactions:
        lines.append("          reactions='all',")
    if haveTransport:
        lines.append("          transport={0!r},".format(transportModel))
    lines.append('          initial_state=state(temperature=300.0, pressure=OneAtm))')
    lines.append('')

    # Write the individual species data
    lines.append(delimiterLine)
    lines.append('# Species data')
    lines.append(delimiterLine)
    lines.append('')

    for s in species:
        lines.append(s.to_cti())

    # Write the reactions
    lines.append(delimiterLine)
    lines.append('# Reaction data')
    lines.append(delimiterLine)

    for i,r in enumerate(reactions):
        lines.append('\n# Reaction {0}'.format(i+1))
        lines.append(r.to_cti())

    lines.append('')

    f = open(outName, 'w')
    f.write('\n'.join(lines))

################################################################################

if __name__ == '__main__':
    import sys
    species, reactions = loadChemkinFile(sys.argv[1])

    if len(sys.argv) > 2:
        lines = open(sys.argv[2]).readlines()
        parseTransportData(lines, species)

    writeCTI(species, reactions)
