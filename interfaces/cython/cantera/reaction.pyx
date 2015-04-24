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
    def __cinit__(self, *args, init=True, **kwargs):
        if init:
            self._reaction.reset(new CxxReaction(0))
            self.reaction = self._reaction.get()

    cdef _assign(self, shared_ptr[CxxReaction] other):
        self._reaction = other
        self.reaction = self._reaction.get()

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

    property products:
        def __get__(self):
            return comp_map_to_dict(self.reaction.products)

    property orders:
        def __get__(self):
            return comp_map_to_dict(self.reaction.orders)

    property ID:
        def __get__(self):
            return pystr(self.reaction.id)

    property reversible:
        def __get__(self):
            return self.reaction.reversible

    property duplicate:
        def __get__(self):
            return self.reaction.duplicate

    property allow_nonreactant_orders:
        def __get__(self):
            return self.reaction.allow_nonreactant_orders

    property allow_negative_orders:
        def __get__(self):
            return self.reaction.allow_negative_orders


cdef class Arrhenius:
    def __cinit__(self, A=0, b=0, E=0, init=True):
        if init:
            self.rate = new CxxArrhenius(A, b, E)
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
                  rate.activationEnergy_R())
    return r


cdef class ElementaryReaction(Reaction):
    property rate:
        def __get__(self):
            cdef CxxElementaryReaction* r = <CxxElementaryReaction*>self.reaction
            return wrapArrhenius(&(r.rate), self)


cdef class ThirdBodyReaction(ElementaryReaction):
    cdef CxxThirdBodyReaction* tbr(self):
        return <CxxThirdBodyReaction*>self.reaction

    property efficiencies:
        def __get__(self):
            return comp_map_to_dict(self.tbr().third_body.efficiencies)

    property default_efficiency:
        def __get__(self):
            return self.tbr().third_body.default_efficiency

    def efficiency(self, species):
        return self.tbr().third_body.efficiency(stringify(species))


cdef class Falloff:
    def __cinit__(self, init=True):
        if init:
            self._falloff.reset(new CxxFalloff())
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


cdef wrapFalloff(shared_ptr[CxxFalloff] falloff):
    f = Falloff(init=False)
    f._falloff = falloff
    f.falloff = f._falloff.get()
    return f


cdef class FalloffReaction(Reaction):
    property low_rate:
        def __get__(self):
            cdef CxxFalloffReaction* r = <CxxFalloffReaction*>self.reaction
            return wrapArrhenius(&(r.low_rate), self)

    property high_rate:
        def __get__(self):
            cdef CxxFalloffReaction* r = <CxxFalloffReaction*>self.reaction
            return wrapArrhenius(&(r.high_rate), self)

    property falloff:
        def __get__(self):
            cdef CxxFalloffReaction* r = <CxxFalloffReaction*>self.reaction
            return wrapFalloff(r.falloff)


cdef class ChemicallyActivatedReaction(FalloffReaction):
    pass


cdef class PlogReaction(Reaction):
    property rates:
        def __get__(self):
            cdef CxxPlogReaction* r = <CxxPlogReaction*>self.reaction
            rates = []
            cdef vector[pair[double,CxxArrhenius]] cxxrates = r.rate.rates()
            cdef pair[double,CxxArrhenius] p_rate
            for p_rate in cxxrates:
                rates.append((p_rate.first,copyArrhenius(&p_rate.second)))
            return rates

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
    else:
        R = Reaction(init=False)

    R._assign(reaction)
    return R
