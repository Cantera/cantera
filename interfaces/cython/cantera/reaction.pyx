cdef extern from "cantera/kinetics/reaction_defs.h" namespace "Cantera":
    cdef int ELEMENTARY_RXN
    cdef int THREE_BODY_RXN
    cdef int FALLOFF_RXN
    cdef int PLOG_RXN
    cdef int CHEBYSHEV_RXN
    cdef int CHEMACT_RXN
    cdef int INTERFACE_RXN

    cdef int SIMPLE_FALLOFF
    cdef int TROE_FALLOFF
    cdef int SRI_FALLOFF


cdef class Reaction:
    reaction_type = 0

    def __cinit__(self, *args, init=True, **kwargs):
        if init:
            self._reaction.reset(newReaction(self.reaction_type))
            self.reaction = self._reaction.get()

    cdef _assign(self, shared_ptr[CxxReaction] other):
        self._reaction = other
        self.reaction = self._reaction.get()

    @staticmethod
    def fromCti(text):
        """
        Create a Reaction object from its CTI string representation.
        """
        cxx_reactions = CxxGetReactions(deref(CxxGetXmlFromString(stringify(text))))
        assert cxx_reactions.size() == 1, cxx_reactions.size()
        return wrapReaction(cxx_reactions[0])

    @staticmethod
    def fromXml(text):
        """
        Create a Reaction object from its XML string representation.
        """
        cxx_reaction = CxxNewReaction(deref(CxxGetXmlFromString(stringify(text))))
        return wrapReaction(cxx_reaction)

    @staticmethod
    def listFromFile(filename):
        """
        Create a list of Reaction objects from all of the reactions defined in a
        CTI or XML file.

        Directories on Cantera's input file path will be searched for the
        specified file.

        In the case of an XML file, the <reactions> nodes are assumed to be
        children of the <reactionsData> node in a document with a <ctml> root
        node, as in the XML files produced by conversion from CTI files.
        """
        cxx_reactions = CxxGetReactions(deref(CxxGetXmlFile(stringify(filename))))
        return [wrapReaction(r) for r in cxx_reactions]

    @staticmethod
    def listFromXml(text):
        """
        Create a list of Reaction objects from all the reaction defined in an
        XML string. The <reaction> nodes are assumed to be children of the
        <reactionData> node in a document with a <ctml> root node, as in the XML
        files produced by conversion from CTI files.
        """
        cxx_reactions = CxxGetReactions(deref(CxxGetXmlFromString(stringify(text))))
        return [wrapReaction(r) for r in cxx_reactions]

    @staticmethod
    def listFromCti(text):
        """
        Create a list of Species objects from all the species defined in a CTI
        string.
        """
        # Currently identical to listFromXml since get_XML_from_string is able
        # to distinguish between CTI and XML.
        cxx_reactions = CxxGetReactions(deref(CxxGetXmlFromString(stringify(text))))
        return [wrapReaction(r) for r in cxx_reactions]

    property reactant_string:
        def __get__(self):
            return pystr(self.reaction.reactantString())

    property product_string:
        def __get__(self):
            return pystr(self.reaction.productString())

    property equation:
        def __get__(self):
            return pystr(self.reaction.equation())

    property reactants:
        def __get__(self):
            return comp_map_to_dict(self.reaction.reactants)
        def __set__(self, reactants):
            self.reaction.reactants = comp_map(reactants)

    property products:
        def __get__(self):
            return comp_map_to_dict(self.reaction.products)
        def __set__(self, products):
            self.reaction.products = comp_map(products)

    property orders:
        def __get__(self):
            return comp_map_to_dict(self.reaction.orders)
        def __set__(self, orders):
            self.reaction.orders = comp_map(orders)

    property ID:
        def __get__(self):
            return pystr(self.reaction.id)
        def __set__(self, ID):
            self.reaction.id = stringify(ID)

    property reversible:
        def __get__(self):
            return self.reaction.reversible
        def __set__(self, reversible):
            self.reaction.reversible = reversible

    property duplicate:
        def __get__(self):
            return self.reaction.duplicate
        def __set__(self, duplicate):
             self.reaction.duplicate = duplicate

    property allow_nonreactant_orders:
        def __get__(self):
            return self.reaction.allow_nonreactant_orders
        def __set__(self, allow):
            self.reaction.allow_nonreactant_orders = allow

    property allow_negative_orders:
        def __get__(self):
            return self.reaction.allow_negative_orders
        def __set__(self, allow):
            self.reaction.allow_negative_orders = allow


cdef class Arrhenius:
    def __cinit__(self, A=0, b=0, E=0, Ta=0, init=True):
        if init:
            self.rate = new CxxArrhenius(A, b, E / gas_constant)
            self.reaction = None

    def __dealloc__(self):
        if self.reaction is None:
            del self.rate

    property preexponential_factor:
        def __get__(self):
            return self.rate.preExponentialFactor()

    property temperature_exponent:
        def __get__(self):
            return self.rate.temperatureExponent()

    property activation_energy:
        def __get__(self):
            return self.rate.activationEnergy_R() * gas_constant


cdef wrapArrhenius(CxxArrhenius* rate, Reaction reaction):
    r = Arrhenius(init=False)
    r.rate = rate
    r.reaction = reaction
    return r

cdef copyArrhenius(CxxArrhenius* rate):
    r = Arrhenius(rate.preExponentialFactor(), rate.temperatureExponent(),
                  rate.activationEnergy_R() * gas_constant)
    return r


cdef class ElementaryReaction(Reaction):
    reaction_type = ELEMENTARY_RXN

    property rate:
        def __get__(self):
            cdef CxxElementaryReaction* r = <CxxElementaryReaction*>self.reaction
            return wrapArrhenius(&(r.rate), self)
        def __set__(self, Arrhenius rate):
            cdef CxxElementaryReaction* r = <CxxElementaryReaction*>self.reaction
            r.rate = deref(rate.rate)


cdef class ThirdBodyReaction(ElementaryReaction):
    reaction_type = THREE_BODY_RXN

    cdef CxxThirdBodyReaction* tbr(self):
        return <CxxThirdBodyReaction*>self.reaction

    property efficiencies:
        def __get__(self):
            return comp_map_to_dict(self.tbr().third_body.efficiencies)
        def __set__(self, eff):
            self.tbr().third_body.efficiencies = comp_map(eff)

    property default_efficiency:
        def __get__(self):
            return self.tbr().third_body.default_efficiency
        def __set__(self, default_eff):
            self.tbr().third_body.default_efficiency = default_eff

    def efficiency(self, species):
        return self.tbr().third_body.efficiency(stringify(species))


cdef class Falloff:
    falloff_type = SIMPLE_FALLOFF

    def __cinit__(self, params=(), init=True):
        if not init:
            return

        cdef vector[double] c
        for p in params:
            c.push_back(p)
        self._falloff = CxxNewFalloff(self.falloff_type, c)
        self.falloff = self._falloff.get()

    property type:
        def __get__(self):
            cdef int falloff_type = self.falloff.getType()
            if falloff_type == SIMPLE_FALLOFF:
                return "Simple"
            elif falloff_type == TROE_FALLOFF:
                return "Troe"
            elif falloff_type == SRI_FALLOFF:
                return "SRI"
            else:
                return "unknown"

    property parameters:
        def __get__(self):
            N = self.falloff.nParameters()
            if N == 0:
                return np.empty(0)

            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(N)
            self.falloff.getParameters(&data[0])
            return data

    def __call__(self, float T, float Pr):
        N = max(self.falloff.workSize(), 1)
        cdef np.ndarray[np.double_t, ndim=1] work = np.empty(N)
        self.falloff.updateTemp(T, &work[0])
        return self.falloff.F(Pr, &work[0])


cdef class TroeFalloff(Falloff):
    falloff_type = TROE_FALLOFF


cdef class SriFalloff(Falloff):
    falloff_type = SRI_FALLOFF


cdef wrapFalloff(shared_ptr[CxxFalloff] falloff):
    cdef int falloff_type = falloff.get().getType()
    if falloff_type == SIMPLE_FALLOFF:
        f = Falloff(init=False)
    elif falloff_type == TROE_FALLOFF:
        f = TroeFalloff(init=False)
    elif falloff_type == SRI_FALLOFF:
        f = SriFalloff(init=False)
    else:
        warnings.warn('Unknown falloff type: {0}'.format(falloff_type))
        f = Falloff(init=False)
    f._falloff = falloff
    f.falloff = f._falloff.get()
    return f


cdef class FalloffReaction(Reaction):
    reaction_type = FALLOFF_RXN

    cdef CxxFalloffReaction* frxn(self):
        return <CxxFalloffReaction*>self.reaction

    property low_rate:
        def __get__(self):
            return wrapArrhenius(&(self.frxn().low_rate), self)
        def __set__(self, Arrhenius rate):
            self.frxn().low_rate = deref(rate.rate)

    property high_rate:
        def __get__(self):
            return wrapArrhenius(&(self.frxn().high_rate), self)
        def __set__(self, Arrhenius rate):
            self.frxn().high_rate = deref(rate.rate)

    property falloff:
        def __get__(self):
            return wrapFalloff(self.frxn().falloff)
        def __set__(self, Falloff f):
            self.frxn().falloff = f._falloff

    property efficiencies:
        def __get__(self):
            return comp_map_to_dict(self.frxn().third_body.efficiencies)
        def __set__(self, eff):
            self.frxn().third_body.efficiencies = comp_map(eff)

    property default_efficiency:
        def __get__(self):
            return self.frxn().third_body.default_efficiency
        def __set__(self, default_eff):
            self.frxn().third_body.default_efficiency = default_eff

    def efficiency(self, species):
        return self.frxn().third_body.efficiency(stringify(species))


cdef class ChemicallyActivatedReaction(FalloffReaction):
    reaction_type = CHEMACT_RXN


cdef class PlogReaction(Reaction):
    reaction_type = PLOG_RXN

    property rates:
        def __get__(self):
            cdef CxxPlogReaction* r = <CxxPlogReaction*>self.reaction
            rates = []
            cdef vector[pair[double,CxxArrhenius]] cxxrates = r.rate.rates()
            cdef pair[double,CxxArrhenius] p_rate
            for p_rate in cxxrates:
                rates.append((p_rate.first,copyArrhenius(&p_rate.second)))
            return rates

        def __set__(self, rates):
            cdef multimap[double,CxxArrhenius] ratemap
            cdef Arrhenius rate
            cdef pair[double,CxxArrhenius] item
            for p,rate in rates:
                item.first = p
                item.second = deref(rate.rate)
                ratemap.insert(item)

            cdef CxxPlogReaction* r = <CxxPlogReaction*>self.reaction
            r.rate = CxxPlog(ratemap)


cdef class ChebyshevReaction(Reaction):
    reaction_type = CHEBYSHEV_RXN

    property Tmin:
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.Tmin()

    property Tmax:
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.Tmax()

    property Pmin:
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.Pmin()

    property Pmax:
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.Pmax()

    property nPressure:
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.nPressure()

    property nTemperature:
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.nTemperature()

    property coeffs:
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            c = np.fromiter(r.rate.coeffs(), np.double)
            return c.reshape((r.rate.nTemperature(), r.rate.nPressure()))


cdef class CoverageDepenency:
    cdef public double a
    cdef public double m
    cdef public double E

    def __init__(self, a, m, E):
        self.a = a
        self.m = m
        self.E = E


cdef class InterfaceReaction(ElementaryReaction):
    reaction_type = INTERFACE_RXN

    property coverage_deps:
        def __get__(self):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            print('ncov:', r.coverage_deps.size())
            deps = {}
            cdef pair[string,CxxCoverageDependency] item
            for item in r.coverage_deps:
                deps[pystr(item.first)] = CoverageDepenency(
                    item.second.a, item.second.m, item.second.E * gas_constant)
            return deps

    property is_sticking_coefficient:
        def __get__(self):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            return r.is_sticking_coefficient

    property sticking_species:
        def __get__(self):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            return pystr(r.sticking_species)


cdef Reaction wrapReaction(shared_ptr[CxxReaction] reaction):
    """
    Wrap a C++ Reaction object with a Python object of the correct derived type.
    """
    cdef int reaction_type = reaction.get().reaction_type

    if reaction_type == ELEMENTARY_RXN:
        R = ElementaryReaction(init=False)
    elif reaction_type == THREE_BODY_RXN:
        R = ThirdBodyReaction(init=False)
    elif reaction_type == FALLOFF_RXN:
        R = FalloffReaction(init=False)
    elif reaction_type == CHEMACT_RXN:
        R = ChemicallyActivatedReaction(init=False)
    elif reaction_type == PLOG_RXN:
        R = PlogReaction(init=False)
    elif reaction_type == CHEBYSHEV_RXN:
        R = ChebyshevReaction(init=False)
    elif reaction_type == INTERFACE_RXN:
        R = InterfaceReaction(init=False)
    else:
        R = Reaction(init=False)

    R._assign(reaction)
    return R

cdef CxxReaction* newReaction(int reaction_type):
    if reaction_type == ELEMENTARY_RXN:
        return new CxxElementaryReaction()
    elif reaction_type == THREE_BODY_RXN:
        return new CxxThirdBodyReaction()
    elif reaction_type == FALLOFF_RXN:
        return new CxxFalloffReaction()
    elif reaction_type == CHEMACT_RXN:
        return new CxxChemicallyActivatedReaction()
    elif reaction_type == PLOG_RXN:
        return new CxxPlogReaction()
    elif reaction_type == CHEBYSHEV_RXN:
        return new CxxChebyshevReaction()
    elif reaction_type == INTERFACE_RXN:
        return new CxxInterfaceReaction()
    else:
        return new CxxReaction(0)
