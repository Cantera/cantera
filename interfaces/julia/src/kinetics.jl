# Kinetics accessors.  All functions accept a `Solution` or a `Kinetics`.

"Number of reactions in the mechanism."
n_reactions(g::KineticsLike) = Int(check(LibCantera.kin_nReactions(_kinetics_handle(g))))

"Total number of species across all phases in the Kinetics manager."
n_total_species(g::KineticsLike) =
    Int(check(LibCantera.kin_nTotalSpecies(_kinetics_handle(g))))

"Reaction equation string for reaction `i` (1-based)."
function reaction_equation(g::KineticsLike, i::Integer)
    r = check(LibCantera.kin_reaction(_kinetics_handle(g), Int32(i - 1)))
    eq = get_string((n, b) -> LibCantera.rxn_equation(r, n, b))
    LibCantera.rxn_del(r)
    return eq
end

"Vector of all reaction equation strings."
reaction_equations(g::KineticsLike) = [reaction_equation(g, i) for i in 1:n_reactions(g)]

"Whether reaction `i` (1-based) is reversible."
is_reversible(g::KineticsLike, i::Integer) =
    check(LibCantera.kin_isReversible(_kinetics_handle(g), Int32(i - 1))) != 0

# Length helpers: rate/ROP arrays are per-reaction; production arrays per-species.
_nrxn(g) = n_reactions(g)
_nsp(g) = n_total_species(g)

# Per-reaction kinetic arrays.
for (jl, c) in (
        (:forward_rates_of_progress, :kin_getFwdRatesOfProgress),
        (:reverse_rates_of_progress, :kin_getRevRatesOfProgress),
        (:net_rates_of_progress,     :kin_getNetRatesOfProgress),
        (:equilibrium_constants,     :kin_getEquilibriumConstants),
        (:forward_rate_constants,    :kin_getFwdRateConstants),
        (:delta_enthalpy,            :kin_getDeltaEnthalpy),
        (:delta_gibbs,               :kin_getDeltaGibbs),
        (:delta_entropy,             :kin_getDeltaEntropy),
    )
    bang = Symbol(jl, :!)
    @eval begin
        $jl(g::KineticsLike) =
            get_array(_nrxn(g), (n, b) -> LibCantera.$c(_kinetics_handle(g), n, b))
        $bang(out, g::KineticsLike) =
            get_array!(out, (n, b) -> LibCantera.$c(_kinetics_handle(g), n, b))
    end
end

"""
    reverse_rate_constants(gas; do_irreversible=false)

Reverse rate constants for all reactions.
"""
function reverse_rate_constants(g::KineticsLike; do_irreversible::Bool=false)
    h = _kinetics_handle(g)
    return get_array(_nrxn(g),
        (n, b) -> LibCantera.kin_getRevRateConstants(h, n, b, Int32(do_irreversible)))
end

# Per-species production arrays.
for (jl, c) in (
        (:net_production_rates, :kin_getNetProductionRates),
        (:creation_rates,       :kin_getCreationRates),
        (:destruction_rates,    :kin_getDestructionRates),
    )
    bang = Symbol(jl, :!)
    @eval begin
        $jl(g::KineticsLike) =
            get_array(_nsp(g), (n, b) -> LibCantera.$c(_kinetics_handle(g), n, b))
        $bang(out, g::KineticsLike) =
            get_array!(out, (n, b) -> LibCantera.$c(_kinetics_handle(g), n, b))
    end
end

# ---- analytic kinetic derivatives -------------------------------------------
# These wrap Cantera's native analytic derivative getters, exposed through the
# CLib via recipes added to interfaces/sourcegen/.../headers/ctkin.yaml (the
# `kin_get*_dd{T,P,C}` functions).  Each is an array getter with the standard
# `(len, ptr)` array convention, so they reuse `get_array`/`get_array!`.
#
# Derivatives are taken at constant volume/composition for `_ddT` and `_ddP`,
# and with respect to molar concentration for `_ddC`, following Cantera's C++
# convention.  The composition Jacobian (`_ddX`) returns a *sparse* matrix in
# C++ and is not representable through the flat-array CLib; it is intentionally
# not wrapped here.  See docs/src/differentiation.md.
#
# Runtime note: these require a `libcantera` built with the ctkin derivative
# recipes (Cantera >= this branch).  Against an older library the underlying
# symbol is absent and the call raises an error at call time.

for (jl, c, len) in (
        (:net_production_rates_ddT, :kin_getNetProductionRates_ddT, :_nsp),
        (:net_production_rates_ddP, :kin_getNetProductionRates_ddP, :_nsp),
        (:net_production_rates_ddC, :kin_getNetProductionRates_ddC, :_nsp),
        (:creation_rates_ddT,       :kin_getCreationRates_ddT,      :_nsp),
        (:creation_rates_ddP,       :kin_getCreationRates_ddP,      :_nsp),
        (:creation_rates_ddC,       :kin_getCreationRates_ddC,      :_nsp),
        (:destruction_rates_ddT,    :kin_getDestructionRates_ddT,   :_nsp),
        (:destruction_rates_ddP,    :kin_getDestructionRates_ddP,   :_nsp),
        (:destruction_rates_ddC,    :kin_getDestructionRates_ddC,   :_nsp),
        (:net_rates_of_progress_ddT, :kin_getNetRatesOfProgress_ddT, :_nrxn),
        (:net_rates_of_progress_ddP, :kin_getNetRatesOfProgress_ddP, :_nrxn),
        (:net_rates_of_progress_ddC, :kin_getNetRatesOfProgress_ddC, :_nrxn),
    )
    bang = Symbol(jl, :!)
    @eval begin
        $jl(g::KineticsLike) =
            get_array($len(g), (n, b) -> LibCantera.$c(_kinetics_handle(g), n, b))
        $bang(out, g::KineticsLike) =
            get_array!(out, (n, b) -> LibCantera.$c(_kinetics_handle(g), n, b))
    end
end

# ---- composition Jacobians (_ddX) -------------------------------------------
# The C++ `_ddX` getters return sparse matrices; the CLib densifies them into a
# flat column-major buffer (see interfaces/sourcegen/.../ctkin.yaml), which we
# reshape into a Julia matrix.

"""
    net_production_rates_ddX(gas) -> Matrix{Float64}

Dense composition Jacobian d(wdot_k)/dX_j of the net production rates with
respect to species mole fractions, size nSpecies x nSpecies.
"""
function net_production_rates_ddX(g::KineticsLike)
    n = _nsp(g)
    buf = get_array(n * n, (len, b) -> LibCantera.kin_getNetProductionRates_ddX(_kinetics_handle(g), len, b))
    return reshape(buf, n, n)
end

"""
    net_rates_of_progress_ddX(gas) -> Matrix{Float64}

Dense composition Jacobian d(q_i)/dX_j of the net rates of progress, size
nReactions x nSpecies.
"""
function net_rates_of_progress_ddX(g::KineticsLike)
    nr = _nrxn(g); ns = _nsp(g)
    buf = get_array(nr * ns, (len, b) -> LibCantera.kin_getNetRatesOfProgress_ddX(_kinetics_handle(g), len, b))
    return reshape(buf, nr, ns)
end

"""
    heat_release_rate(gas) -> Float64

Volumetric heat release rate [W/m^3], `-Σ_k h_k · ẇ_k`, where `h_k` are the
partial molar enthalpies and `ẇ_k` the net production rates.  Matches Python's
`gas.heat_release_rate` (the CLib exposes no direct getter).
"""
heat_release_rate(g::Solution) =
    -sum(partial_molar_enthalpies(g) .* net_production_rates(g))
