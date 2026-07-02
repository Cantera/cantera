# The Solution type: the primary user-facing object, mirroring Python's
# `cantera.Solution`.  A Solution bundles a ThermoPhase, a Kinetics manager and
# a Transport manager over a single phase, and owns their CLib handles.

"""
    Solution(infile; name="", transport="")
    Solution(infile, name; transport="")

Construct a `Solution` from a Cantera input (YAML) file.

`name` selects the phase within the file (the first phase is used when empty);
`transport` selects a transport model (the file's default when empty).

```julia
gas = Solution("gri30.yaml")
gas = Solution("gri30.yaml", "gri30")
```
"""
mutable struct Solution <: CanteraObject
    handle::Int32
    thermo::Int32
    kinetics::Int32
    transport::Int32
    closed::Bool

    function Solution(infile::AbstractString, name::AbstractString="";
                      transport::AbstractString="default")
        h = check(LibCantera.sol_newSolution(infile, name, transport))
        th = check(LibCantera.sol_thermo(h))
        # Kinetics/Transport may be absent for some phases; tolerate failures.
        kin = _try_handle(() -> LibCantera.sol_kinetics(h))
        tr = _try_handle(() -> LibCantera.sol_transport(h))
        s = new(h, th, kin, tr, false)
        finalizer(close!, s)
        return s
    end
end

"Run a handle-returning CLib call, returning -1 instead of throwing on failure."
function _try_handle(f)
    h = try
        f()
    catch
        return Int32(-1)
    end
    return (h == -1 || h == _ERR) ? Int32(-1) : Int32(h)
end

"""
    close!(gas::Solution)

Release the Cantera-side objects backing `gas`.  Safe to call multiple times;
the finalizer calls this automatically, so explicit use is only needed for
deterministic cleanup.
"""
function close!(s::Solution)
    s.closed && return s
    s.closed = true
    # Best-effort deletion; never throw from a finalizer path.
    for (h, del) in ((s.transport, LibCantera.trans_del),
                     (s.kinetics, LibCantera.kin_del),
                     (s.thermo, LibCantera.thermo_del),
                     (s.handle, LibCantera.sol_del))
        h == -1 && continue
        try
            del(h)
        catch
        end
    end
    return s
end

function Base.show(io::IO, s::Solution)
    if s.closed
        print(io, "Solution(<closed>)")
    else
        print(io, "Solution(", n_species(s), " species, ",
              n_reactions(s), " reactions)")
    end
end

# ---- handle accessors -------------------------------------------------------
# Thermo/kinetics/transport functions accept either a Solution or the specific
# sub-phase wrapper.  These helpers extract the right handle for dispatch.

_thermo_handle(s::Solution) = s.thermo
_thermo_handle(t::ThermoPhase) = t.handle

function _kinetics_handle(s::Solution)
    s.kinetics == -1 &&
        throw(CanteraError("this Solution has no Kinetics manager"))
    return s.kinetics
end
_kinetics_handle(k::Kinetics) = k.handle

function _transport_handle(s::Solution)
    s.transport == -1 &&
        throw(CanteraError("this Solution has no Transport manager"))
    return s.transport
end
_transport_handle(t::Transport) = t.handle

const ThermoLike = Union{Solution,ThermoPhase}
const KineticsLike = Union{Solution,Kinetics}
const TransportLike = Union{Solution,Transport}

"""
    thermo(gas::Solution) -> ThermoPhase

Return the [`ThermoPhase`](@ref) view of `gas`.
"""
thermo(s::Solution) = ThermoPhase(s.thermo)

"""
    kinetics(gas::Solution) -> Kinetics

Return the [`Kinetics`](@ref) view of `gas`.
"""
kinetics(s::Solution) = Kinetics(_kinetics_handle(s))

"""
    transport(gas::Solution) -> Transport

Return the [`Transport`](@ref) view of `gas`.
"""
transport(s::Solution) = Transport(_transport_handle(s))

"""
    name(gas::Solution) -> String

Name of the Solution object.
"""
name(s::Solution) = get_string((n, b) -> LibCantera.sol_name(s.handle, n, b))

"""
    transport_model(gas::Solution) -> String

Name of the active transport model.
"""
transport_model(s::Solution) =
    get_string((n, b) -> LibCantera.sol_transportModel(s.handle, n, b))
