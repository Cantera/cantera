# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# One-dimensional flames (Domain1D / Sim1D).
#
# Parity wrapper over the `ctdomain` and `ctonedim` CLib libraries, mirroring
# Python/MATLAB `FreeFlame`.  A [`FreeFlame`](@ref) owns three [`Domain1D`](@ref)
# objects (inlet / flow / outlet) plus a `Sim1D` solver handle, and keeps a
# reference to the gas [`Solution`](@ref) alive for their lifetime.

# ---------------------------------------------------------------------------
# Domain1D: a generic wrapper around a single 1-D domain handle.
# ---------------------------------------------------------------------------

"""
    Domain1D

Wrapper around a single Cantera 1-D domain (a boundary such as an inlet/outlet,
or a flow region).  Domains are normally created and managed by a
[`FreeFlame`](@ref); the accessors here let you inspect and set profiles on an
individual domain.  Domains borrow the gas [`Solution`](@ref) handle and are
freed together with the owning simulation.
"""
mutable struct Domain1D <: CanteraObject
    handle::Int32
    solution::Solution   # keep the phase alive
    closed::Bool
end

"Number of solution points (grid nodes) in the domain."
n_points(d::Domain1D) = Int(check(LibCantera.domain_nPoints(d.handle)))

"Number of solution components in the domain."
n_components(d::Domain1D) = Int(check(LibCantera.domain_nComponents(d.handle)))

"Name of component `n` (1-based)."
function component_name(d::Domain1D, n::Integer)
    return get_string(
        (bl, b) -> LibCantera.domain_componentName(d.handle, Int32(n - 1), bl, b))
end

"Vector of all component names in the domain."
component_names(d::Domain1D) = [component_name(d, n) for n in 1:n_components(d)]

"Domain type string (e.g. \"free-flow\", \"inlet\", \"outlet\")."
domain_type(d::Domain1D) =
    get_string((bl, b) -> LibCantera.domain_domainType(d.handle, bl, b))

"Grid points of the domain [m]."
grid(d::Domain1D) =
    get_array(n_points(d), (bl, b) -> LibCantera.domain_grid(d.handle, bl, b))

"Scalar value of `component` in a boundary domain (a single-point domain)."
value(d::Domain1D, component::AbstractString) =
    checkd(LibCantera.domain_value(d.handle, component))

"Per-point profile of `component` across the domain grid."
function solution_profile(d::Domain1D, component::AbstractString)
    np = n_points(d)
    return get_array(np,
        (n, b) -> LibCantera.domain_getValues(d.handle, component, n, b))
end

"""
    set_profile!(domain, component, positions, values)

Set the profile of `component` over normalized positions `positions` (in
`[0, 1]`, spanning the domain) with the given `values`.
"""
function set_profile!(d::Domain1D, component::AbstractString, positions, values)
    pos = as_f64(collect(positions))
    val = as_f64(collect(values))
    check(LibCantera.domain_setProfile(d.handle, component,
                                       Int32(length(pos)), pointer(pos),
                                       Int32(length(val)), pointer(val)))
    return d
end

"Set a spatially uniform value of `component` across the whole domain."
function set_flat_profile!(d::Domain1D, component::AbstractString, value)
    check(LibCantera.domain_setFlatProfile(d.handle, component, Float64(value)))
    return d
end

"Set a uniform grid of `points` nodes spanning `[start, start+length]` [m]."
function setup_uniform_grid!(d::Domain1D, points::Integer, length::Real,
                              start::Real=0.0)
    check(LibCantera.domain_setupUniformGrid(d.handle, Int32(points), Float64(length),
                                             Float64(start)))
    return d
end

"Set an explicit grid from a vector of positions [m]."
function setup_grid!(d::Domain1D, z)
    zz = as_f64(collect(z))
    check(LibCantera.domain_setupGrid(d.handle, Int32(length(zz)), pointer(zz)))
    return d
end

function Base.show(io::IO, d::Domain1D)
    if d.closed
        print(io, "Domain1D(<closed>)")
    else
        print(io, "Domain1D(\"", domain_type(d), "\", ", n_points(d), " points)")
    end
end

# ---------------------------------------------------------------------------
# FreeFlame: a freely-propagating premixed flame.
# ---------------------------------------------------------------------------

"""
    FreeFlame(gas; width=0.03, points=8)

Assemble a freely-propagating premixed flame from the current state of `gas`,
mirroring Python's `cantera.FreeFlame`.  The unburned mixture state (temperature,
pressure, composition) is taken from `gas` at construction time.

The flame consists of three domains in solver order: an inlet boundary, a
`free-flow` region of the given `width` [m] discretized with `points` initial
grid nodes, and an outlet boundary.  Call [`solve!`](@ref) to compute the
solution; the laminar flame speed is then available from [`flame_speed`](@ref).
"""
mutable struct FreeFlame <: CanteraObject
    gas::Solution
    inlet::Domain1D
    flow::Domain1D
    outlet::Domain1D
    sim::Int32
    rho_u::Float64          # unburned density [kg/m^3], for the flame-speed eigenvalue
    Tu::Float64             # unburned temperature [K]
    P::Float64              # pressure [Pa]
    Yu::Vector{Float64}     # unburned mass fractions
    closed::Bool
end

# Initial (non-uniform) flame grid used by Python's FreeFlame for a given width.
_initial_flame_grid(width) = as_f64([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0] .* width)

function FreeFlame(gas::Solution; width::Real=0.03)
    # Capture the unburned state from the incoming gas.
    Tu    = temperature(gas)
    P     = pressure(gas)
    Xu    = mole_fractions(gas)
    Yu    = mass_fractions(gas)
    rho_u = density(gas)

    # Create the three domains (order matters: inlet, flow, outlet).
    hin  = check(LibCantera.domain_newBoundary1D("inlet",  gas.handle, "reactants"))
    hflw = check(LibCantera.domain_newFlow1D("free-flow",  gas.handle, "flame"))
    hout = check(LibCantera.domain_newBoundary1D("outlet", gas.handle, "products"))

    inlet  = Domain1D(hin,  gas, false)
    flow   = Domain1D(hflw, gas, false)
    outlet = Domain1D(hout, gas, false)

    # Flow domain: operating pressure and the initial non-uniform flame grid.
    check(LibCantera.flow_setPressure(hflw, Float64(P)))
    z0 = _initial_flame_grid(width)
    GC.@preserve z0 check(LibCantera.domain_setupGrid(hflw, Int32(length(z0)),
                                                       pointer(z0)))

    # Inlet boundary: unburned temperature and composition.  The mass flux is
    # left at zero here; `set_initial_guess!` seeds it (mdot = 1*rho_u) exactly
    # as Python does.
    check(LibCantera.bdry_setTemperature(hin, Float64(Tu)))
    xu = as_f64(Xu)
    GC.@preserve xu check(LibCantera.bdry_setMoleFractions(hin, Int32(length(xu)),
                                                            pointer(xu)))

    # Build the simulation from the three domain handles (order matters).
    # GC.@preserve keeps `doms` alive across the ccall.s
    doms = Int32[hin, hflw, hout]
    hsim = GC.@preserve doms check(LibCantera.sim1D_newSim1D(Int32(length(doms)),
                                                              pointer(doms)))

    flame = FreeFlame(gas, inlet, flow, outlet, hsim, rho_u, Tu, P, Yu, false)
    finalizer(close!, flame)

    _set_initial_guess!(flame)
    return flame
end

# Replicates Python `FreeFlame.set_initial_guess` (locs = [0, 0.3, 0.5, 1.0]).
function _set_initial_guess!(flame::FreeFlame)
    hflw = flame.flow.handle
    gas  = flame.gas
    locs = as_f64([0.0, 0.3, 0.5, 1.0])

    _setprof(comp, vals) = GC.@preserve locs vals check(
        LibCantera.domain_setProfile(hflw, comp, Int32(length(locs)), pointer(locs),
                                     Int32(length(vals)), pointer(vals)))

    # Unburned state at the inlet.
    set_TPY!(gas, flame.Tu, flame.P, flame.Yu)

    mdot = checkd(LibCantera.bdry_mdot(flame.inlet.handle))
    if mdot == 0.0
        # A nonzero initial guess increases the likelihood of convergence.
        mdot = 1.0 * density(gas)
        check(LibCantera.bdry_setMdot(flame.inlet.handle, mdot))
    end

    Y0 = flame.Yu
    u0 = mdot / density(gas)
    T0 = flame.Tu

    # Burned (adiabatic-flame) state by equilibrating at constant H, P.
    equilibrate!(gas, "HP")
    Teq = temperature(gas)
    Yeq = mass_fractions(gas)
    u1  = mdot / density(gas)

    _setprof("velocity", as_f64([u0, u0, u1, u1]))
    _setprof("T", as_f64([T0, T0, Teq, Teq]))

    # Pick the fixed-temperature point, reusing an existing grid point when a
    # reasonable choice exists (Python's exact criterion).
    T = solution_profile(flame.flow, "T")
    Tmid = 0.75 * T0 + 0.25 * Teq
    i = findlast(<(Tmid), T)          # 1-based analogue of numpy flatnonzero[-1]
    Tfixed = if Tmid - T[i] < 0.2 * (Tmid - T0)
        T[i]
    elseif T[i + 1] - Tmid < 0.2 * (Teq - Tmid)
        T[i + 1]
    else
        Tmid
    end
    check(LibCantera.sim1D_setFixedTemperature(flame.sim, Float64(Tfixed)))

    # Species mass-fraction profiles: unburned -> burned.
    for n in 1:n_species(gas)
        _setprof(species_name(gas, n), as_f64([Y0[n], Y0[n], Yeq[n], Yeq[n]]))
    end
    return flame
end

"""
    set_refine_criteria!(flame; ratio=10.0, slope=0.8, curve=0.8, prune=0.0)

Set the grid-refinement criteria on the flow domain (domain index 1).  Defaults
match Python's `FlameBase.set_refine_criteria`.
"""
function set_refine_criteria!(flame::FreeFlame; ratio=10.0, slope=0.8, curve=0.8,
                               prune=0.0)
    check(LibCantera.sim1D_setRefineCriteria(flame.sim, Int32(1), Float64(ratio),
                                             Float64(slope), Float64(curve),
                                             Float64(prune)))
    return flame
end

"""
    set_inlet!(flame; T=nothing, X=nothing, mdot=nothing)

Update the inlet (reactants) boundary temperature, composition and/or mass flux.
"""
function set_inlet!(flame::FreeFlame; T=nothing, X=nothing, mdot=nothing)
    h = flame.inlet.handle
    T === nothing || check(LibCantera.bdry_setTemperature(h, Float64(T)))
    mdot === nothing || check(LibCantera.bdry_setMdot(h, Float64(mdot)))
    if X !== nothing
        if X isa AbstractString
            check(LibCantera.bdry_setMoleFractionsByName(h, X))
        else
            x = as_f64(collect(X))
            check(LibCantera.bdry_setMoleFractions(h, Int32(length(x)), pointer(x)))
        end
    end
    return flame
end

"Fix the flame temperature at `T` [K] (used by the free-flame eigenvalue solve)."
function set_fixed_temperature!(flame::FreeFlame, T::Real)
    check(LibCantera.sim1D_setFixedTemperature(flame.sim, Float64(T)))
    return flame
end

# Regrid the flow domain to a uniform `N`-point grid over `[zmin, zmax]`.
function _regrid!(flame::FreeFlame, N::Integer, zmin::Real, zmax::Real)
    z = collect(range(zmin, zmax; length=N))
    GC.@preserve z check(LibCantera.domain_setupGrid(flame.flow.handle,
                                                     Int32(length(z)), pointer(z)))
    return flame
end

# One `sim1D_solve`; returns `true` on success, `false` on a CanteraError.
function _try_solve(sim, ll, refine)
    try
        check(LibCantera.sim1D_solve(sim, ll, Int32(refine)))
        return true
    catch
        return false
    end
end

# Signals that the flow domain is too narrow to contain the flame; caught by
# `solve!`, which doubles the domain width (mirrors Python's `DomainTooNarrow`).
struct _DomainTooNarrow <: Exception end

# Python's `check_width` criterion: the domain is too narrow if the temperature
# gradient at either edge is significant compared to the average gradient.
function _width_ok(flame::FreeFlame)
    T = solution_profile(flame.flow, "T")
    x = grid(flame.flow)
    length(x) < 3 && return true
    mRef = (T[end] - T[1]) / (x[end] - x[1])
    mRef == 0 && return true
    mLeft = (T[2] - T[1]) / (x[2] - x[1]) / mRef
    mRight = (T[end-2] - T[end]) / (x[end-2] - x[end]) / mRef
    return !(mLeft > 0.02 || mRight > 0.02)
end

# Replicates the Cython `Sim1D.solve(auto=True)` staged loop: progressively
# finer uniform grids (8, 12, 24, 48 points), re-seeding the initial guess on
# each, trying energy-on first and falling back to a fixed-temperature solve,
# then a refinement pass.  Returns `true` once a solution is found.
function _staged_solve!(flame::FreeFlame, ll::Int32, refine_grid::Bool)
    sim  = flame.sim
    flow = flame.flow.handle
    g = grid(flame.flow)
    zmin, zmax = g[1], g[end]

    nPoints = [n_points(flame.flow), 12, 24, 48]
    solved = false
    for N in nPoints
        if N != n_points(flame.flow)
            _regrid!(flame, N, zmin, zmax)
        end
        _set_initial_guess!(flame)

        # Try with the energy equation enabled first (usually works).
        check(LibCantera.flow_setEnergyEnabled(flow, Int32(1)))
        solved = _try_solve(sim, ll, false)

        # Fall back to a fixed-temperature solve, then re-enable energy.
        if !solved
            check(LibCantera.flow_setEnergyEnabled(flow, Int32(0)))
            solved = _try_solve(sim, ll, false)
            if solved
                check(LibCantera.flow_setEnergyEnabled(flow, Int32(1)))
                solved = _try_solve(sim, ll, false)
            end
        end

        # Emulate Python's steady-state width callback: it fires after this
        # (pre-refinement) steady solve on the coarse grid, so check here.
        solved && !_width_ok(flame) && throw(_DomainTooNarrow())

        # Refinement pass; a non-extinct converged solution ends the loop.
        # (`extinct()` is always false for a free flame.)
        if solved && refine_grid
            solved = _try_solve(sim, ll, true)
            solved && !_width_ok(flame) && throw(_DomainTooNarrow())
            solved && break
        end

        refine_grid || break
    end
    solved || throw(CanteraError("Could not find a solution for the 1D problem"))
    return solved
end

"""
    solve!(flame; loglevel=0, refine_grid=true, auto=true)

Solve the flame.

With `auto=true` (the default) this reproduces Python's
`FreeFlame.solve(auto=True)` exactly: an internal staged multi-grid schedule
wrapped in a domain-widening loop.  After each staged solve the temperature
gradients at the domain edges are checked; if the flame is too close to a
boundary the grid is doubled (and refined) and the staged solve is repeated,
up to 12 times.

!!! note
    Python additionally installs the width check as a *steady-state callback*
    inside the C++ solver, so it can abort and widen mid-solve.  The CLib does
    not expose steady callbacks, so here the identical width criterion is applied
    *after* each completed staged solve instead.  This is the one piece of the
    reference algorithm that is emulated externally rather than in-solver; the
    final converged result is unaffected for the usual case.

With `auto=false` a single `sim1D_solve(loglevel, refine_grid)` is issued using
the current refine criteria (advanced use).
"""
function solve!(flame::FreeFlame; loglevel::Integer=0, refine_grid::Bool=true,
                 auto::Bool=true)
    flame.closed && throw(CanteraError("FreeFlame is closed"))
    ll = Int32(loglevel)
    sim = flame.sim
    flow = flame.flow.handle

    if !auto
        if !_try_solve(sim, ll, refine_grid)
            throw(CanteraError("Could not find a solution for the 1D problem"))
        end
        return flame
    end

    # Domain-widening loop (mirrors FreeFlame.solve check_width / grid *= 2).
    for _ in 1:12
        try
            _staged_solve!(flame, ll, refine_grid)
            break
        catch err
            err isa _DomainTooNarrow || rethrow()
            # Too narrow: double the domain, then re-run the *staged* solve from
            # scratch. Mirrors Python's `self.flame.grid *= 2; self.refine(...)`:
            # `sim1D_refine` only adapts the grid (inserting points per the
            # refine criteria), it does not solve.
            znew = grid(flame.flow) .* 2
            GC.@preserve znew check(LibCantera.domain_setupGrid(flow,
                                        Int32(length(znew)), pointer(znew)))
            refine_grid && check(LibCantera.sim1D_refine(sim, ll))
        end
    end
    return flame
end

# ---- post-solve accessors ---------------------------------------------------

"Grid points of the flow (flame) domain [m]."
grid(flame::FreeFlame) = grid(flame.flow)

"Per-point profile of `component` in the flow domain."
solution_profile(flame::FreeFlame, component::AbstractString) =
    solution_profile(flame.flow, component)

"Temperature profile across the flame [K]."
flame_T(flame::FreeFlame) = solution_profile(flame.flow, "T")

"Axial velocity profile across the flame [m/s]."
flame_velocity(flame::FreeFlame) = solution_profile(flame.flow, "velocity")

"Mole-fraction profile of `species` across the flame."
flame_X(flame::FreeFlame, species::AbstractString) =
    solution_profile(flame.flow, species)

"""
    flame_speed(flame) -> Float64

Laminar flame speed [m/s]: the inlet (unburned) axial velocity.  For a freely
propagating flame the inlet mass flux is an eigenvalue updated by the solver, so
this reads the current inlet `mdot` and divides by the unburned density.
"""
flame_speed(flame::FreeFlame) =
    checkd(LibCantera.bdry_mdot(flame.inlet.handle)) / flame.rho_u

"Number of grid points in the flow domain."
n_points(flame::FreeFlame) = n_points(flame.flow)

function close!(flame::FreeFlame)
    flame.closed && return flame
    flame.closed = true
    try
        LibCantera.sim1D_del(flame.sim)
    catch
    end
    for d in (flame.inlet, flame.flow, flame.outlet)
        d.closed = true
        try
            LibCantera.domain_del(d.handle)
        catch
        end
    end
    return flame
end

function Base.show(io::IO, flame::FreeFlame)
    if flame.closed
        print(io, "FreeFlame(<closed>)")
    else
        print(io, "FreeFlame(", n_points(flame), " grid points)")
    end
end

# ---------------------------------------------------------------------------
# BurnerFlame: a burner-stabilized flat flame.
# ---------------------------------------------------------------------------

"""
    BurnerFlame(gas; width=0.03, mdot=nothing)

Assemble a burner-stabilized flat flame from the current state of `gas`,
mirroring Python's `cantera.BurnerFlame`.  The unburned mixture state
(temperature, pressure, composition) is taken from `gas` at construction time.

Unlike a [`FreeFlame`](@ref), the burner inlet has a user-prescribed mass flux
`mdot` [kg/m^2/s]; there is no flame-speed eigenvalue and no fixed-temperature
anchor.  If `mdot` is not given it defaults to `0.4 * ρ_unburned`.

The flame consists of three domains in solver order: a burner (inlet) boundary,
an `unstrained-flow` region of the given `width` [m], and an outlet boundary.
Call [`solve!`](@ref) to compute the solution.  Use [`set_burner!`](@ref) to
change the burner temperature, composition, or mass flux.
"""
mutable struct BurnerFlame <: CanteraObject
    gas::Solution
    burner::Domain1D
    flow::Domain1D
    outlet::Domain1D
    sim::Int32
    rho_u::Float64          # unburned density [kg/m^3]
    Tu::Float64             # unburned temperature [K]
    P::Float64              # pressure [Pa]
    Yu::Vector{Float64}     # unburned mass fractions
    closed::Bool
end

# Initial (non-uniform) grid used by Python's BurnerFlame for a given width.
_initial_burner_grid(width) = as_f64([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0] .* width)

function BurnerFlame(gas::Solution; width::Real=0.03, mdot=nothing)
    # Capture the unburned (burner) state from the incoming gas.
    Tu    = temperature(gas)
    P     = pressure(gas)
    Xu    = mole_fractions(gas)
    Yu    = mass_fractions(gas)
    rho_u = density(gas)

    # Create the three domains (order matters: burner, flow, outlet).
    hbrn = check(LibCantera.domain_newBoundary1D("inlet", gas.handle, "burner"))
    hflw = check(LibCantera.domain_newFlow1D("unstrained-flow", gas.handle, "flame"))
    hout = check(LibCantera.domain_newBoundary1D("outlet", gas.handle, "products"))

    burner = Domain1D(hbrn, gas, false)
    flow   = Domain1D(hflw, gas, false)
    outlet = Domain1D(hout, gas, false)

    # Flow domain: operating pressure and the initial non-uniform grid.
    check(LibCantera.flow_setPressure(hflw, Float64(P)))
    z0 = _initial_burner_grid(width)
    GC.@preserve z0 check(LibCantera.domain_setupGrid(hflw, Int32(length(z0)),
                                                       pointer(z0)))

    # Burner boundary: unburned temperature, composition, and prescribed mdot.
    check(LibCantera.bdry_setTemperature(hbrn, Float64(Tu)))
    xu = as_f64(Xu)
    GC.@preserve xu check(LibCantera.bdry_setMoleFractions(hbrn, Int32(length(xu)),
                                                            pointer(xu)))
    md = mdot === nothing ? 0.4 * rho_u : Float64(mdot)
    check(LibCantera.bdry_setMdot(hbrn, md))

    doms = Int32[hbrn, hflw, hout]
    hsim = GC.@preserve doms check(LibCantera.sim1D_newSim1D(Int32(length(doms)),
                                                              pointer(doms)))

    flame = BurnerFlame(gas, burner, flow, outlet, hsim, rho_u, Tu, P, Yu, false)
    finalizer(close!, flame)

    _set_initial_guess!(flame)
    return flame
end

# Replicates Python `BurnerFlame.set_initial_guess` (locs = [0, 0.2, 1.0]).
function _set_initial_guess!(flame::BurnerFlame)
    hflw = flame.flow.handle
    gas  = flame.gas
    locs = as_f64([0.0, 0.2, 1.0])

    _setprof(comp, vals) = GC.@preserve locs vals check(
        LibCantera.domain_setProfile(hflw, comp, Int32(length(locs)), pointer(locs),
                                     Int32(length(vals)), pointer(vals)))

    # Unburned state at the burner.
    set_TPY!(gas, flame.Tu, flame.P, flame.Yu)
    mdot = checkd(LibCantera.bdry_mdot(flame.burner.handle))

    Y0 = flame.Yu
    u0 = mdot / density(gas)
    T0 = flame.Tu

    # Burned (adiabatic-flame) state by equilibrating at constant H, P.
    equilibrate!(gas, "HP")
    Teq = temperature(gas)
    Yeq = mass_fractions(gas)
    u1  = mdot / density(gas)

    _setprof("velocity", as_f64([u0, u1, u1]))
    _setprof("T", as_f64([T0, Teq, Teq]))
    for n in 1:n_species(gas)
        _setprof(species_name(gas, n), as_f64([Y0[n], Yeq[n], Yeq[n]]))
    end
    return flame
end

"""
    set_burner!(flame; T=nothing, X=nothing, mdot=nothing)

Update the burner (reactants) boundary temperature, composition and/or mass flux
[kg/m^2/s].
"""
function set_burner!(flame::BurnerFlame; T=nothing, X=nothing, mdot=nothing)
    h = flame.burner.handle
    if T !== nothing
        check(LibCantera.bdry_setTemperature(h, Float64(T)))
        flame.Tu = Float64(T)
    end
    mdot === nothing || check(LibCantera.bdry_setMdot(h, Float64(mdot)))
    if X !== nothing
        if X isa AbstractString
            check(LibCantera.bdry_setMoleFractionsByName(h, X))
        else
            x = as_f64(collect(X))
            check(LibCantera.bdry_setMoleFractions(h, Int32(length(x)), pointer(x)))
        end
        # Refresh the cached unburned composition/density from the new state.
        set_TPX!(flame.gas, flame.Tu, flame.P, X isa AbstractString ? X : collect(X))
        flame.Yu = mass_fractions(flame.gas)
        flame.rho_u = density(flame.gas)
    end
    return flame
end

"""
    set_refine_criteria!(flame::BurnerFlame; ratio=10.0, slope=0.8, curve=0.8,
                          prune=0.0)

Set the grid-refinement criteria on the flow domain (domain index 1).
"""
function set_refine_criteria!(flame::BurnerFlame; ratio=10.0, slope=0.8, curve=0.8,
                               prune=0.0)
    check(LibCantera.sim1D_setRefineCriteria(flame.sim, Int32(1), Float64(ratio),
                                             Float64(slope), Float64(curve),
                                             Float64(prune)))
    return flame
end

_regrid!(flame::BurnerFlame, N::Integer, zmin::Real, zmax::Real) =
    (z = collect(range(zmin, zmax; length=N));
     GC.@preserve z check(LibCantera.domain_setupGrid(flame.flow.handle,
                                                      Int32(length(z)), pointer(z)));
     flame)

# Staged multi-grid solve for a burner flame: no fixed-temperature eigenvalue.
function _staged_solve!(flame::BurnerFlame, ll::Int32, refine_grid::Bool)
    sim  = flame.sim
    flow = flame.flow.handle
    g = grid(flame.flow)
    zmin, zmax = g[1], g[end]

    nPoints = [n_points(flame.flow), 12, 24, 48]
    solved = false
    for N in nPoints
        if N != n_points(flame.flow)
            _regrid!(flame, N, zmin, zmax)
        end
        _set_initial_guess!(flame)

        check(LibCantera.flow_setEnergyEnabled(flow, Int32(1)))
        solved = _try_solve(sim, ll, false)
        if !solved
            check(LibCantera.flow_setEnergyEnabled(flow, Int32(0)))
            solved = _try_solve(sim, ll, false)
            if solved
                check(LibCantera.flow_setEnergyEnabled(flow, Int32(1)))
                solved = _try_solve(sim, ll, false)
            end
        end
        if solved && refine_grid
            solved = _try_solve(sim, ll, true)
            solved && break
        end
        refine_grid || break
    end
    solved || throw(CanteraError("Could not find a solution for the 1D problem"))
    return solved
end

"""
    solve!(flame::BurnerFlame; loglevel=0, refine_grid=true, auto=true)

Solve the burner-stabilized flame.  With `auto=true` a staged multi-grid
schedule is used; with `auto=false` a single `sim1D_solve` is issued using the
current refine criteria.  Unlike a free flame there is no flame-speed eigenvalue
or fixed-temperature anchor; the burner mass flux is a fixed input.
"""
function solve!(flame::BurnerFlame; loglevel::Integer=0, refine_grid::Bool=true,
                 auto::Bool=true)
    flame.closed && throw(CanteraError("BurnerFlame is closed"))
    ll = Int32(loglevel)
    if !auto
        if !_try_solve(flame.sim, ll, refine_grid)
            throw(CanteraError("Could not find a solution for the 1D problem"))
        end
        return flame
    end
    _staged_solve!(flame, ll, refine_grid)
    return flame
end

# ---- post-solve accessors ---------------------------------------------------

grid(flame::BurnerFlame) = grid(flame.flow)
solution_profile(flame::BurnerFlame, component::AbstractString) =
    solution_profile(flame.flow, component)
flame_T(flame::BurnerFlame) = solution_profile(flame.flow, "T")
flame_velocity(flame::BurnerFlame) = solution_profile(flame.flow, "velocity")
flame_X(flame::BurnerFlame, species::AbstractString) =
    solution_profile(flame.flow, species)
n_points(flame::BurnerFlame) = n_points(flame.flow)

"Prescribed burner mass flux [kg/m^2/s]."
burner_mdot(flame::BurnerFlame) = checkd(LibCantera.bdry_mdot(flame.burner.handle))

function close!(flame::BurnerFlame)
    flame.closed && return flame
    flame.closed = true
    try
        LibCantera.sim1D_del(flame.sim)
    catch
    end
    for d in (flame.burner, flame.flow, flame.outlet)
        d.closed = true
        try
            LibCantera.domain_del(d.handle)
        catch
        end
    end
    return flame
end

function Base.show(io::IO, flame::BurnerFlame)
    if flame.closed
        print(io, "BurnerFlame(<closed>)")
    else
        print(io, "BurnerFlame(", n_points(flame), " grid points)")
    end
end

export Domain1D, FreeFlame, BurnerFlame, set_burner!, burner_mdot,
       n_points, n_components, component_name,
       component_names, domain_type, grid, value, solution_profile, set_profile!,
       set_flat_profile!, setup_uniform_grid!, setup_grid!, set_refine_criteria!,
       solve!, flame_speed, flame_T, flame_X, flame_velocity,
       set_fixed_temperature!, set_inlet!
