# Reactor wrappers.  A Reactor holds a handle into the CLib reactor cabinet and
# a reference to the Solution it was built from (to keep it alive).

"""
    Reactor

Wrapper around a Cantera reactor CLib handle.  Construct one of the concrete
helpers ([`IdealGasReactor`](@ref), [`Reactor`](@ref), [`ConstPressureReactor`](@ref),
[`IdealGasConstPressureReactor`](@ref)) rather than this type directly.
"""
mutable struct Reactor <: CanteraObject
    handle::Int32
    solution::Solution   # keep the phase alive for the reactor's lifetime
    closed::Bool
end

function _new_reactor(model::AbstractString, gas::Solution, name::AbstractString)
    h = check(LibCantera.reactor_new(model, gas.handle, Int32(0), name))
    r = Reactor(h, gas, false)
    finalizer(close!, r)
    return r
end

"""
    IdealGasReactor(gas; name="")

Zero-dimensional constant-volume reactor with an ideal-gas energy equation.
"""
IdealGasReactor(gas::Solution; name::AbstractString="") =
    _new_reactor("IdealGasReactor", gas, name)

"""
    Reactor(gas; name="")

General constant-volume reactor.
"""
Reactor(gas::Solution; name::AbstractString="") =
    _new_reactor("Reactor", gas, name)

"""
    ConstPressureReactor(gas; name="")

Constant-pressure reactor.
"""
ConstPressureReactor(gas::Solution; name::AbstractString="") =
    _new_reactor("ConstPressureReactor", gas, name)

"""
    IdealGasConstPressureReactor(gas; name="")

Constant-pressure reactor with an ideal-gas energy equation.
"""
IdealGasConstPressureReactor(gas::Solution; name::AbstractString="") =
    _new_reactor("IdealGasConstPressureReactor", gas, name)

function close!(r::Reactor)
    r.closed && return r
    r.closed = true
    try
        LibCantera.reactor_del(r.handle)
    catch
    end
    return r
end

# ---- scalar reactor state ---------------------------------------------------

for (jl, c) in (
        (:mass,                :reactor_mass),
        (:volume,              :reactor_volume),
        (:density,             :reactor_density),
        (:temperature,         :reactor_temperature),
        (:pressure,            :reactor_pressure),
        (:enthalpy_mass,       :reactor_enthalpy_mass),
        (:internal_energy_mass,:reactor_intEnergy_mass),
    )
    @eval $jl(r::Reactor) = checkd(LibCantera.$c(r.handle))
end

"Reactor name."
name(r::Reactor) = get_string((n, b) -> LibCantera.reactor_name(r.handle, n, b))

"Reactor type string."
reactor_type(r::Reactor) = get_string((n, b) -> LibCantera.reactor_type(r.handle, n, b))

"Mass fractions in the reactor."
function mass_fractions(r::Reactor)
    nsp = n_species(r.solution)
    return get_array(nsp, (n, b) -> LibCantera.reactor_massFractions(r.handle, n, b))
end

"Enable/disable chemistry in the reactor."
function set_energy_enabled!(r::Reactor, flag::Bool)
    check(LibCantera.reactor_setEnergyEnabled(r.handle, Int32(flag)))
    return r
end

"Enable/disable the energy equation in the reactor."
function set_chemistry_enabled!(r::Reactor, flag::Bool)
    check(LibCantera.reactor_setChemistryEnabled(r.handle, Int32(flag)))
    return r
end

"Set the initial reactor volume [m^3]."
function set_initial_volume!(r::Reactor, vol)
    check(LibCantera.reactor_setInitialVolume(r.handle, Float64(vol)))
    return r
end

Base.show(io::IO, r::Reactor) =
    print(io, r.closed ? "Reactor(<closed>)" : "$(reactor_type(r))(\"$(name(r))\")")

# ---- reservoirs -------------------------------------------------------------

"""
    Reservoir(gas; name="")

A reactor with a fixed thermodynamic state, used as a boundary (source or sink)
for [`FlowDevice`](@ref)s and [`Wall`](@ref)s.  Its state never changes during
integration.
"""
Reservoir(gas::Solution; name::AbstractString="") =
    _new_reactor("Reservoir", gas, name)

# ---- surfaces ---------------------------------------------------------------

"""
    reactor_area(r) -> Float64

Wall/surface area associated with the reactor [m^2].
"""
area(r::Reactor) = checkd(LibCantera.reactor_area(r.handle))

"Set the wall/surface area of the reactor [m^2]."
function set_area!(r::Reactor, a)
    check(LibCantera.reactor_setArea(r.handle, Float64(a)))
    return r
end

"""
    ReactorSurface(surf, reactor; name="")

A surface on which heterogeneous reactions take place, coupling the surface
phase `surf` (a [`Solution`](@ref) for an interface phase) to a bulk `reactor`.
"""
function ReactorSurface(surf::Solution, reactor::Reactor; name::AbstractString="")
    reactors = Int32[reactor.handle]
    h = check(LibCantera.reactor_newSurface(surf.handle, Int32(length(reactors)),
                                            pointer(reactors), Int32(0), name))
    # Reuse the Reactor wrapper: a ReactorSurface is a ReactorBase in the CLib
    # cabinet and is deleted with reactor_del.  Keep both phases alive.
    rs = Reactor(h, surf, false)
    finalizer(close!, rs)
    return rs
end

# ---- mass flow rate (for flow reactors / surfaces) --------------------------

"Mass flow rate through the reactor [kg/s]."
mass_flow_rate(r::Reactor) = checkd(LibCantera.reactor_massFlowRate(r.handle))

"Set the mass flow rate through the reactor [kg/s]."
function set_mass_flow_rate!(r::Reactor, mdot)
    check(LibCantera.reactor_setMassFlowRate(r.handle, Float64(mdot)))
    return r
end

# ---- sensitivities ----------------------------------------------------------

"""
    add_sensitivity_reaction!(r, i)

Mark reaction `i` (1-based) as a sensitivity parameter for reactor `r`.  The
reactor must already be part of a [`ReactorNet`](@ref).
"""
function add_sensitivity_reaction!(r::Reactor, i::Integer)
    check(LibCantera.reactor_addSensitivityReaction(r.handle, Int32(i - 1)))
    return r
end

"Number of sensitivity parameters associated with the reactor."
n_sens_params(r::Reactor) = Int(check(LibCantera.reactor_nSensParams(r.handle)))
