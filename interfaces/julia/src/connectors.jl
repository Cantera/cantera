# Reactor connectors: walls and flow devices.
#
# A connector links two reactors in a network.  Each Julia wrapper stores the
# CLib connector handle and keeps references to both reactors alive for its
# lifetime.  Connectors are freed with `connector_del` via a finalizer.

"""
    Connector

Abstract supertype of objects that link two reactors in a network
([`Wall`](@ref), [`MassFlowController`](@ref), [`Valve`](@ref),
[`PressureController`](@ref)).
"""
abstract type Connector <: CanteraObject end

"""
    FlowDevice <: Connector

Abstract supertype of flow devices that move mass between reactors.
"""
abstract type FlowDevice <: Connector end

# Resolve a Func1-like object (or raw handle) to its CLib handle.
_func_handle(f::Integer) = Int32(f)
_func_handle(f) = Int32(f.handle)

# ---- concrete types ---------------------------------------------------------

for T in (:Wall, :MassFlowController, :Valve, :PressureController)
    super = T === :Wall ? :Connector : :FlowDevice
    @eval begin
        mutable struct $T <: $super
            handle::Int32
            r0::Reactor      # upstream / left
            r1::Reactor      # downstream / right
            closed::Bool
        end
    end
end

function _new_connector(::Type{T}, model, left::Reactor, right::Reactor,
                        name) where {T<:Connector}
    h = check(LibCantera.connector_new(model, left.handle, right.handle, name))
    c = T(h, left, right, false)
    finalizer(close!, c)
    return c
end

function close!(c::Connector)
    c.closed && return c
    c.closed = true
    try
        LibCantera.connector_del(c.handle)
    catch
    end
    return c
end

"Connector name."
name(c::Connector) = get_string((n, b) -> LibCantera.connector_name(c.handle, n, b))

"Connector type string."
connector_type(c::Connector) =
    get_string((n, b) -> LibCantera.connector_type(c.handle, n, b))

"Set the connector name."
function set_name!(c::Connector, nm::AbstractString)
    check(LibCantera.connector_setName(c.handle, nm))
    return c
end

# ---- walls ------------------------------------------------------------------

"""
    Wall(left, right; A=1.0, U=0.0, K=0.0, Q=nothing, velocity=nothing,
         expansion_rate_coeff=nothing, emissivity=nothing, name="")

A wall separating reactors `left` and `right`.  `A` is the wall area [m^2], `U`
the heat-transfer coefficient [W/m^2/K], `K` the expansion-rate coefficient
[m/s/Pa].  `Q` (heat flux) and `velocity` may be `Func1` objects; `emissivity`
enables radiative transfer.
"""
function Wall(left::Reactor, right::Reactor; A=1.0, U=0.0, K=0.0, Q=nothing,
              velocity=nothing, expansion_rate_coeff=nothing, emissivity=nothing,
              name::AbstractString="")
    w = _new_connector(Wall, "Wall", left, right, name)
    set_area!(w, A)
    U != 0 && set_heat_transfer_coeff!(w, U)
    K != 0 && set_expansion_rate_coeff!(w, K)
    expansion_rate_coeff !== nothing && set_expansion_rate_coeff!(w, expansion_rate_coeff)
    emissivity !== nothing && set_emissivity!(w, emissivity)
    Q !== nothing && set_heat_flux!(w, Q)
    velocity !== nothing && set_velocity!(w, velocity)
    return w
end

"Rate of volumetric expansion of the left reactor [m^3/s]."
expansion_rate(w::Wall) = checkd(LibCantera.wall_expansionRate(w.handle))

"Rate of heat transfer through the wall (left to right) [W]."
heat_rate(w::Wall) = checkd(LibCantera.wall_heatRate(w.handle))

"Wall area [m^2]."
area(w::Wall) = checkd(LibCantera.wall_area(w.handle))

"Set the wall area [m^2]."
function set_area!(w::Wall, a)
    check(LibCantera.wall_setArea(w.handle, Float64(a)))
    return w
end

"Set the wall heat-transfer coefficient [W/m^2/K]."
function set_heat_transfer_coeff!(w::Wall, U)
    check(LibCantera.wall_setHeatTransferCoeff(w.handle, Float64(U)))
    return w
end

"Set the wall thermal resistance [m^2*K/W]."
function set_thermal_resistance!(w::Wall, Rth)
    check(LibCantera.wall_setThermalResistance(w.handle, Float64(Rth)))
    return w
end

"Set the wall expansion-rate coefficient [m/s/Pa]."
function set_expansion_rate_coeff!(w::Wall, k)
    check(LibCantera.wall_setExpansionRateCoeff(w.handle, Float64(k)))
    return w
end

"Set the wall emissivity for radiative heat transfer (0..1)."
function set_emissivity!(w::Wall, epsilon)
    check(LibCantera.wall_setEmissivity(w.handle, Float64(epsilon)))
    return w
end

"Set the wall heat flux as a time-dependent `Func1` (or handle) [W/m^2]."
function set_heat_flux!(w::Wall, q)
    check(LibCantera.wall_setHeatFlux(w.handle, _func_handle(q)))
    return w
end

"Set the wall velocity as a time-dependent `Func1` (or handle) [m/s]."
function set_velocity!(w::Wall, f)
    check(LibCantera.wall_setVelocity(w.handle, _func_handle(f)))
    return w
end

# ---- flow devices -----------------------------------------------------------

"Mass flow rate through the flow device [kg/s]."
mass_flow_rate(d::FlowDevice) = checkd(LibCantera.flowdev_massFlowRate(d.handle))

"Device coefficient of the flow device (meaning depends on the device type)."
device_coefficient(d::FlowDevice) = checkd(LibCantera.flowdev_deviceCoefficient(d.handle))

"Set the device coefficient of the flow device."
function set_device_coefficient!(d::FlowDevice, c)
    check(LibCantera.flowdev_setDeviceCoefficient(d.handle, Float64(c)))
    return d
end

"Set the pressure function of the flow device as a `Func1` (or handle)."
function set_pressure_function!(d::FlowDevice, f)
    check(LibCantera.flowdev_setPressureFunction(d.handle, _func_handle(f)))
    return d
end

"Set the time function of the flow device as a `Func1` (or handle)."
function set_time_function!(d::FlowDevice, g)
    check(LibCantera.flowdev_setTimeFunction(d.handle, _func_handle(g)))
    return d
end

"Set the primary flow device (used by [`PressureController`](@ref))."
function set_primary!(d::FlowDevice, primary::FlowDevice)
    check(LibCantera.flowdev_setPrimary(d.handle, primary.handle))
    return d
end

"""
    MassFlowController(upstream, downstream; mdot=0.0, name="")

A flow device that maintains a specified mass flow rate `mdot` [kg/s] from
`upstream` to `downstream`, independent of the pressure difference.
"""
function MassFlowController(upstream::Reactor, downstream::Reactor; mdot=0.0,
                            name::AbstractString="")
    d = _new_connector(MassFlowController, "MassFlowController", upstream, downstream, name)
    mdot != 0 && set_mass_flow_rate!(d, mdot)
    return d
end

"""
    set_mass_flow_rate!(mfc, mdot)

Set the (constant) mass flow rate of a [`MassFlowController`](@ref) [kg/s].
"""
set_mass_flow_rate!(d::MassFlowController, mdot) = set_device_coefficient!(d, mdot)

"""
    Valve(upstream, downstream; K=0.0, name="")

A flow device whose mass flow rate is proportional to the pressure difference,
`mdot = K * (P_upstream - P_downstream)`, with valve coefficient `K`.
"""
function Valve(upstream::Reactor, downstream::Reactor; K=0.0, name::AbstractString="")
    d = _new_connector(Valve, "Valve", upstream, downstream, name)
    K != 0 && set_device_coefficient!(d, K)
    return d
end

"""
    PressureController(upstream, downstream; primary=nothing, K=0.0, name="")

A flow device that regulates the pressure difference across it relative to a
`primary` flow device: `mdot = primary.mdot + K * (P_upstream - P_downstream)`.
"""
function PressureController(upstream::Reactor, downstream::Reactor; primary=nothing,
                           K=0.0, name::AbstractString="")
    d = _new_connector(PressureController, "PressureController", upstream, downstream, name)
    K != 0 && set_device_coefficient!(d, K)
    primary !== nothing && set_primary!(d, primary)
    return d
end

# ---- display ----------------------------------------------------------------

for T in (:Wall, :MassFlowController, :Valve, :PressureController)
    @eval Base.show(io::IO, c::$T) =
        print(io, c.closed ? string($(QuoteNode(T)), "(<closed>)") :
              string($(QuoteNode(T)), "(\"", name(c), "\")"))
end
