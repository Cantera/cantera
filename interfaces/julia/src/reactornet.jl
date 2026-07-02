# ReactorNet: time integration of a network of reactors.

"""
    ReactorNet(reactors::Vector{Reactor})
    ReactorNet(reactor::Reactor)

Create a reactor network for time integration.  The network keeps references to
its reactors so they are not finalized while integration is in progress.
"""
mutable struct ReactorNet <: CanteraObject
    handle::Int32
    reactors::Vector{Reactor}
    closed::Bool

    function ReactorNet(reactors::Vector{Reactor})
        isempty(reactors) && throw(ArgumentError("ReactorNet needs at least one reactor"))
        handles = Int32[r.handle for r in reactors]
        h = check(LibCantera.reactornet_new(Int32(length(handles)), pointer(handles)))
        net = new(h, reactors, false)
        finalizer(close!, net)
        return net
    end
end

ReactorNet(r::Reactor) = ReactorNet([r])

function close!(net::ReactorNet)
    net.closed && return net
    net.closed = true
    try
        LibCantera.reactornet_del(net.handle)
    catch
    end
    return net
end

"""
    advance!(net, t)

Advance the network state to absolute time `t` [s].
"""
function advance!(net::ReactorNet, t)
    check(LibCantera.reactornet_advance(net.handle, Float64(t)))
    return net
end

"""
    step!(net) -> Float64

Take one internal timestep and return the new time [s].
"""
step!(net::ReactorNet) = checkd(LibCantera.reactornet_step(net.handle))

"Current network time [s].  Extends `Base.time` for `ReactorNet`."
Base.time(net::ReactorNet) = checkd(LibCantera.reactornet_time(net.handle))

"Set the initial integration time [s]."
function set_initial_time!(net::ReactorNet, t)
    check(LibCantera.reactornet_setInitialTime(net.handle, Float64(t)))
    return net
end

"Set the maximum internal timestep [s]."
function set_max_time_step!(net::ReactorNet, dt)
    check(LibCantera.reactornet_setMaxTimeStep(net.handle, Float64(dt)))
    return net
end

"""
    set_tolerances!(net; rtol=1e-9, atol=1e-15)

Set the relative and absolute integrator tolerances.
"""
function set_tolerances!(net::ReactorNet; rtol=1e-9, atol=1e-15)
    check(LibCantera.reactornet_setTolerances(net.handle, Float64(rtol), Float64(atol)))
    return net
end

"Relative error tolerance of the network integrator."
rtol(net::ReactorNet) = checkd(LibCantera.reactornet_rtol(net.handle))

"Absolute error tolerance of the network integrator."
atol(net::ReactorNet) = checkd(LibCantera.reactornet_atol(net.handle))

"""
    set_sensitivity_tolerances!(net; rtol=1e-6, atol=1e-6)

Set the relative and absolute tolerances used for sensitivity analysis.
"""
function set_sensitivity_tolerances!(net::ReactorNet; rtol=1e-6, atol=1e-6)
    check(LibCantera.reactornet_setSensitivityTolerances(net.handle,
                                                         Float64(rtol), Float64(atol)))
    return net
end

"""
    sensitivity(net, component, p, reactor) -> Float64

Normalized sensitivity of `component` (e.g. `"temperature"` or a species name)
in `reactor` with respect to sensitivity parameter `p` (1-based).  `reactor` may
be a [`Reactor`](@ref) belonging to the network or its 1-based position.
"""
function sensitivity(net::ReactorNet, component::AbstractString, p::Integer,
                     reactor::Reactor)
    idx = findfirst(r -> r === reactor, net.reactors)
    idx === nothing && throw(ArgumentError("reactor is not part of this network"))
    return sensitivity(net, component, p, idx)
end

function sensitivity(net::ReactorNet, component::AbstractString, p::Integer,
                     reactor::Integer)
    return checkd(LibCantera.reactornet_sensitivity(net.handle, component,
                                                    Int32(p - 1), Int32(reactor - 1)))
end

"""
    n_components(net) -> Int

Number of state variables (equations) integrated by the network.
"""
n_components(net::ReactorNet) = Int(check(LibCantera.reactornet_neq(net.handle)))

"""
    state(net) -> Vector{Float64}

Current network state vector, of length [`n_components`](@ref).
"""
state(net::ReactorNet) =
    get_array(n_components(net), (n, b) -> LibCantera.reactornet_getState(net.handle, n, b))

"Name of state-vector component `i` (1-based)."
function component_name(net::ReactorNet, i::Integer)
    return get_string((n, b) -> LibCantera.reactornet_componentName(net.handle, Int32(i - 1), n, b))
end

"""
    component_names(net) -> Vector{String}

Names of all state-vector components, aligned with [`state`](@ref).
"""
component_names(net::ReactorNet) = [component_name(net, i) for i in 1:n_components(net)]

Base.show(io::IO, net::ReactorNet) =
    print(io, net.closed ? "ReactorNet(<closed>)" :
          "ReactorNet($(length(net.reactors)) reactor(s), t=$(time(net)) s)")
