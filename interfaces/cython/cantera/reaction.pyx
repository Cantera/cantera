cdef extern from "cantera/kinetics/reaction_defs.h" namespace "Cantera":
    cdef int ELEMENTARY_RXN
    cdef int THREE_BODY_RXN
    cdef int FALLOFF_RXN
    cdef int PLOG_RXN
    cdef int CHEBYSHEV_RXN
    cdef int CHEMACT_RXN
    cdef int INTERFACE_RXN



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
    def __cinit__(self, init=True):
        if init:
            self.rate = new CxxArrhenius()
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


cdef Reaction wrapReaction(shared_ptr[CxxReaction] reaction):
    """
    Wrap a C++ Reaction object with a Python object of the correct derived type.
    """
    cdef int reaction_type = reaction.get().reaction_type

    if reaction_type == ELEMENTARY_RXN:
        R = ElementaryReaction(init=False)
    elif reaction_type == THREE_BODY_RXN:
        R = ThirdBodyReaction(init=False)
    else:
        R = Reaction(init=False)

    R._assign(reaction)
    return R
