# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Error handling for the Cantera CLib boundary.
#
# CLib functions signal failure with sentinel return values:
#   * `int32_t` functions return -1 (or ERR = -999);
#   * `double`  functions return DERR = -999.999.
# When such a sentinel is seen, the human-readable message is retrieved with
# `ct_getCanteraError` and raised as a `CanteraError`.

const _ERR = Int32(-999)     # ERR   from clib_defs.h
const _DERR = -999.999       # DERR  from clib_defs.h

"""
    CanteraError(msg)

Exception raised when an underlying Cantera CLib call reports an error.  The
message is the exception text produced by the C++ Cantera core.
"""
struct CanteraError <: Exception
    msg::String
end

Base.showerror(io::IO, e::CanteraError) = print(io, "CanteraError: ", e.msg)

export CanteraError

"Retrieve and clear the last error message stored by the Cantera CLib."
function last_cantera_error()
    buflen = LibCantera.ct_getCanteraError(Int32(0), Ptr{UInt8}(C_NULL))
    buflen <= 0 && return "unknown Cantera error (no message available)"
    buf = Vector{UInt8}(undef, buflen)
    LibCantera.ct_getCanteraError(Int32(buflen), pointer(buf))
    return unsafe_string(pointer(buf))
end

"""
    check(code::Integer) -> Int32

Verify the return code of an `int32_t` CLib call and raise [`CanteraError`](@ref)
on the error sentinel.  Returns the code unchanged on success so it can be used
inline, e.g. `n = check(LibCantera.thermo_nSpecies(h))`.
"""
function check(code::Integer)
    # All successful `int32_t` CLib returns are non-negative (handles, counts,
    # indices, string lengths, or the 0 status).  Any negative value signals an
    # exception.  Callers that must distinguish a "not found" npos result (also
    # negative) do so *before* calling `check` (see `species_index`).
    code < 0 && throw(CanteraError(last_cantera_error()))
    return Int32(code)
end

"""
    checkd(value::Float64) -> Float64

Verify the return value of a `double` CLib call and raise [`CanteraError`](@ref)
on the `DERR` sentinel.
"""
function checkd(value::Float64)
    value == _DERR && throw(CanteraError(last_cantera_error()))
    return value
end
