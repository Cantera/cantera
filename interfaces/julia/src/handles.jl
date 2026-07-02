# Wrapper types and low-level handle helpers.
#
# Every Cantera CLib object is referenced by an `Int32` handle into a C++-side
# "cabinet".  The Julia wrapper types below store that handle and are the only
# things user code touches; raw handles are never part of the public API.

"""
    CanteraObject

Abstract supertype of all Julia wrappers around Cantera CLib objects.
"""
abstract type CanteraObject end

# Sub-phase wrappers.  These *borrow* handles owned by a [`Solution`](@ref); they
# do not carry finalizers of their own (the owning Solution frees everything).
"""
    ThermoPhase

Thermodynamic-phase view of a [`Solution`](@ref).  Wraps a `ThermoPhase` CLib
handle.
"""
mutable struct ThermoPhase <: CanteraObject
    handle::Int32
end

"""
    Kinetics

Kinetics view of a [`Solution`](@ref).  Wraps a `Kinetics` CLib handle.
"""
mutable struct Kinetics <: CanteraObject
    handle::Int32
end

"""
    Transport

Transport view of a [`Solution`](@ref).  Wraps a `Transport` CLib handle.
"""
mutable struct Transport <: CanteraObject
    handle::Int32
end

"`handle(obj)` returns the raw CLib handle backing a wrapper (internal)."
handle(o::CanteraObject) = o.handle

# ---- generic string / array marshalling ------------------------------------

"""
    get_string(f) -> String

Call a CLib string getter through the closure `f(bufLen::Int32, buf::Ptr{UInt8})`
using the size-query protocol: first with a null buffer to obtain the required
length, then again with a correctly sized buffer.
"""
function get_string(f)
    n = check(f(Int32(0), Ptr{UInt8}(C_NULL)))
    n <= 0 && return ""
    buf = Vector{UInt8}(undef, n)
    check(f(Int32(n), pointer(buf)))
    return unsafe_string(pointer(buf))
end

"""
    get_array!(out::Vector{Float64}, f) -> out

Fill `out` in place through the closure `f(len::Int32, ptr::Ptr{Float64})` and
return it.  Used by the performance-oriented `foo!(out, gas)` variants.
"""
function get_array!(out::Vector{Float64}, f)
    check(f(Int32(length(out)), pointer(out)))
    return out
end

"`get_array(n, f)` allocates a `Vector{Float64}` of length `n` and fills it."
get_array(n::Integer, f) = get_array!(Vector{Float64}(undef, n), f)

"Convert an arbitrary real vector to a contiguous `Vector{Float64}` at the boundary."
as_f64(x::Vector{Float64}) = x
as_f64(x::AbstractVector{<:Real}) = convert(Vector{Float64}, x)
