# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

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

"""
    ThermoPhase

Thermodynamic-phase view of a [`Solution`](@ref).  Wraps a `ThermoPhase` CLib
handle and holds a reference to the owning `Solution`.
"""
mutable struct ThermoPhase{P} <: CanteraObject
    handle::Int32
    parent::P
end

"""
    Kinetics

Kinetics view of a [`Solution`](@ref).  Wraps a `Kinetics` CLib handle and holds
a reference to the owning `Solution`.
"""
mutable struct Kinetics{P} <: CanteraObject
    handle::Int32
    parent::P
end

"""
    Transport

Transport view of a [`Solution`](@ref).  Wraps a `Transport` CLib handle and
holds a reference to the owning `Solution`.
"""
mutable struct Transport{P} <: CanteraObject
    handle::Int32
    parent::P
end

"`handle(obj)` returns the raw CLib handle backing a wrapper (internal)."
handle(o::CanteraObject) = o.handle

export CanteraObject, ThermoPhase, Kinetics, Transport

# ---- generic string / array marshalling ------------------------------------

"""
    get_string(f) -> String

Call a CLib string getter through the closure `f(bufLen::Int32, buf)`
using the size-query protocol: first with a null buffer to obtain the required
length, then again with a correctly sized buffer.
"""
function get_string(f)
    n = check(f(Int32(0), Ptr{UInt8}(C_NULL)))
    n <= 0 && return ""
    buf = Vector{UInt8}(undef, n)
    check(f(Int32(n), buf))
    return GC.@preserve buf unsafe_string(pointer(buf))
end

"Whether the elements of `a` are stored contiguously in column-major order."
is_contiguous(a::Array) = true
function is_contiguous(a::StridedArray)
    expected = 1
    for (st, sz) in zip(strides(a), size(a))
        sz == 1 || st == expected || return false
        expected *= sz
    end
    return true
end

"""
    get_array!(out, n, f) -> out

Fill `out` in place through the closure `f(len::Int32, buf)`, which hands `buf`
(`out` itself) to the CLib, and return it.  Used by the performance-oriented
`foo!(out, gas)` variants.

`out` must hold exactly `n` `Float64` values in contiguous storage: a `Vector`,
a `Matrix`, or a unit-stride view such as `view(M, :, j)`.  Throws a
`DimensionMismatch` for any other length, rather than leaving trailing elements
untouched, and an `ArgumentError` for a non-contiguous view.
"""
function get_array!(out::StridedArray{Float64}, n::Integer, f)
    length(out) == n || throw(DimensionMismatch(
        "output array has length $(length(out)), expected $n"))
    is_contiguous(out) || throw(ArgumentError(
        "output array must be contiguous in memory, e.g. `view(M, :, j)`"))
    check(f(Int32(n), out))
    return out
end

"""
    get_matrix!(out, nrows, ncols, f) -> out

Matrix counterpart of [`get_array!`](@ref): fill the `nrows × ncols` matrix `out`
from a CLib getter that writes a flat column-major buffer of that size.
"""
function get_matrix!(out::StridedMatrix{Float64}, nrows::Integer, ncols::Integer, f)
    size(out) == (nrows, ncols) || throw(DimensionMismatch(
        "output matrix has size $(size(out)), expected ($nrows, $ncols)"))
    return get_array!(out, nrows * ncols, f)
end

"""
    inplace_doc(name, arg="gas") -> String

Docstring of the in-place variant `name!(out, arg)` of the allocating getter
`name`.
"""
inplace_doc(name, arg="gas") = """
        $(name)!(out, $arg) -> out

    In-place [`$name`](@ref): fill `out` and return it, without allocating.
    `out` must have the same size as the result of `$name`; a contiguous view
    such as `view(M, :, j)` is accepted.
    """

"`get_array(n, f)` allocates a `Vector{Float64}` of length `n` and fills it."
get_array(n::Integer, f) = get_array!(Vector{Float64}(undef, n), n, f)

"Convert an arbitrary real vector to a contiguous `Vector{Float64}` at the boundary."
as_f64(x::Vector{Float64}) = x
as_f64(x::AbstractVector{<:Real}) = convert(Vector{Float64}, x)
