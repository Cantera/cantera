# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Func1 -- functors of a single variable.
#
# Cantera's `Func1` objects wrap C++ function objects (sin, cos, polynomials,
# compound sum/product functions, ...) that can be evaluated, differentiated
# symbolically and combined arithmetically.  Each object is referenced by an
# `Int32` handle into the `func1` cabinet.

"""
    Func1

A wrapper around a Cantera `Func1` function object.  A `Func1` is
callable: `f(t)` evaluates it at `t`.  Instances combine with `+`, `-`, `*` and
`/` to build compound functions and support symbolic [`derivative`](@ref).

# Examples
```julia
f = Func1("sin", 2.0)   # t -> sin(2t)
f(pi/4)                 # == sin(pi/2)
g = Func1("polynomial3", [1.0, 2.0, 3.0, 4.0])  # t^3 + 2t^2 + 3t + 4
h = Func1("constant", 2.0) + Func1("constant", 3.0)
```
"""
mutable struct Func1 <: CanteraObject
    handle::Int32
    closed::Bool

    function Func1(handle::Integer)
        f = new(Int32(handle), false)
        finalizer(close!, f)
        return f
    end
end

"""
    Func1(type::String, coeff::Real=1.0)

Construct a basic functor of the given `type` (e.g. `"sin"`, `"cos"`, `"exp"`,
`"log"`, `"pow"`, `"constant"`).  For `"sin"` with coefficient `w`, evaluating
at `t` returns `sin(w*t)`.
"""
Func1(type::AbstractString, coeff::Real=1.0) =
    Func1(check(LibCantera.func1_newBasic(type, Float64(coeff))))

"""
    Func1(type::String, coeffs::AbstractVector)

Construct an advanced functor parametrized by an array of coefficients, e.g.
`"polynomial3"` (coefficients from highest to lowest degree) or `"Fourier"`.
"""
function Func1(type::AbstractString, coeffs::AbstractVector{<:Real})
    c = as_f64(coeffs)
    return Func1(check(
        LibCantera.func1_newAdvanced(type, Int32(length(c)), pointer(c))))
end

"""
    constant_function(c) -> Func1

Convenience constructor for the constant functor `t -> c`.
"""
constant_function(c::Real) = Func1("constant", c)

"""
    close!(f::Func1)

Release the Cantera-side object backing `f`.  Idempotent; called automatically
by the finalizer.
"""
function close!(f::Func1)
    f.closed && return f
    f.closed = true
    try
        LibCantera.func1_del(f.handle)
    catch
    end
    return f
end

"""
    evaluate(f::Func1, t) -> Float64

Evaluate the functor at `t`.
"""
evaluate(f::Func1, t::Real) = checkd(LibCantera.func1_eval(f.handle, Float64(t)))

(f::Func1)(t::Real) = evaluate(f, t)

"""
    derivative(f::Func1) -> Func1

Return a new `Func1` representing the symbolic derivative of `f`.
"""
derivative(f::Func1) = Func1(check(LibCantera.func1_derivative(f.handle)))

"""
    func_type(f::Func1) -> String

The functor's type string (e.g. `"sin"`, `"sum"`).
"""
func_type(f::Func1) = get_string((n, b) -> LibCantera.func1_type(f.handle, n, b))

"""
    write(f::Func1, arg="t") -> String

Return a string (LaTeX-style) representation of `f` using `arg` as the
independent-variable symbol.
"""
write(f::Func1, arg::AbstractString="t") =
    get_string((n, b) -> LibCantera.func1_write(f.handle, arg, n, b))

# Arithmetic combinations build compound functors.
Base.:+(f::Func1, g::Func1) =
    Func1(check(LibCantera.func1_newSumFunction(f.handle, g.handle)))
Base.:-(f::Func1, g::Func1) =
    Func1(check(LibCantera.func1_newDiffFunction(f.handle, g.handle)))
Base.:*(f::Func1, g::Func1) =
    Func1(check(LibCantera.func1_newProdFunction(f.handle, g.handle)))
Base.:/(f::Func1, g::Func1) =
    Func1(check(LibCantera.func1_newRatioFunction(f.handle, g.handle)))

function Base.show(io::IO, f::Func1)
    if f.closed
        print(io, "Func1(<closed>)")
    else
        print(io, "Func1(\"", func_type(f), "\": ", write(f), ")")
    end
end

export Func1, evaluate, derivative, func_type, constant_function
