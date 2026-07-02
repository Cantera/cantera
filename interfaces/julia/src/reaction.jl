# Reaction object wrapper.  Reactions are usually accessed through a Kinetics
# manager (see `reaction_equations`), but a standalone wrapper is provided for
# completeness and parity with the CLib `ctrxn` library.

"""
    Reaction

Wrapper around a Cantera `Reaction` CLib handle.
"""
mutable struct Reaction <: CanteraObject
    handle::Int32
    closed::Bool

    function Reaction(handle::Integer)
        r = new(Int32(handle), false)
        finalizer(close!, r)
        return r
    end
end

"""
    reaction(gas, i) -> Reaction

Return reaction `i` (1-based) from a Kinetics manager as a [`Reaction`](@ref).
"""
function reaction(g::KineticsLike, i::Integer)
    h = check(LibCantera.kin_reaction(_kinetics_handle(g), Int32(i - 1)))
    return Reaction(h)
end

function close!(r::Reaction)
    r.closed && return r
    r.closed = true
    try
        LibCantera.rxn_del(r.handle)
    catch
    end
    return r
end

"Reaction equation string."
equation(r::Reaction) = get_string((n, b) -> LibCantera.rxn_equation(r.handle, n, b))

"Reaction type string (e.g. `\"reaction\"`, `\"three-body\"`)."
reaction_type(r::Reaction) = get_string((n, b) -> LibCantera.rxn_type(r.handle, n, b))

"Whether the reaction involves a third body."
uses_third_body(r::Reaction) = check(LibCantera.rxn_usesThirdBody(r.handle)) != 0

Base.show(io::IO, r::Reaction) =
    print(io, "Reaction(\"", r.closed ? "<closed>" : equation(r), "\")")
