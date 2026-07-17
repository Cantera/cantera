# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# MultiPhase mixtures.
#
# A `MultiPhase` holds one or more phases (each a `Solution`) together with the
# number of moles of each, and computes multi-phase equilibrium and mixture
# thermodynamic properties.  It mirrors Python's `cantera.Mixture`.

"""
    MultiPhase()
    MultiPhase(phases::AbstractVector{<:Pair{Solution,<:Real}})

A mixture of one or more phases.  Construct empty and add phases with
[`add_phase!`](@ref) followed by [`init!`](@ref), or pass a vector of
`solution => moles` pairs to add and initialize in one step.

```julia
gas = Solution("gri30.yaml")
mix = MultiPhase([gas => 1.0])
set_temperature!(mix, 1200.0); set_pressure!(mix, one_atm)
equilibrate!(mix, "TP")
```
"""
mutable struct MultiPhase <: CanteraObject
    handle::Int32
    phases::Vector{Solution}   # keep added Solutions alive
    initialized::Bool
    closed::Bool

    function MultiPhase()
        h = check(LibCantera.mix_new())
        mp = new(h, Solution[], false, false)
        finalizer(close!, mp)
        return mp
    end
end

function MultiPhase(phases::AbstractVector{<:Pair{Solution,<:Real}})
    mp = MultiPhase()
    for (sol, moles) in phases
        add_phase!(mp, sol, moles)
    end
    init!(mp)
    return mp
end

"""
    close!(mp::MultiPhase)

Release the Cantera-side object backing `mp`.  Idempotent.
"""
function close!(mp::MultiPhase)
    mp.closed && return mp
    mp.closed = true
    try
        LibCantera.mix_del(mp.handle)
    catch
    end
    return mp
end

"""
    add_phase!(mp::MultiPhase, solution::Solution, moles::Real)

Add `solution` to the mixture with the given number of `moles`.  Call
[`init!`](@ref) after all phases have been added.
"""
function add_phase!(mp::MultiPhase, sol::Solution, moles::Real)
    check(LibCantera.mix_addPhase(mp.handle, sol.handle, Float64(moles)))
    push!(mp.phases, sol)
    return mp
end

"""
    init!(mp::MultiPhase)

Finalize the mixture after all phases have been added.  Required before any
property or equilibrium calculation.
"""
function init!(mp::MultiPhase)
    check(LibCantera.mix_init(mp.handle))
    mp.initialized = true
    return mp
end

"""
    update_phases!(mp::MultiPhase)

Synchronize the constituent phase objects with the mixture's current state.
"""
function update_phases!(mp::MultiPhase)
    check(LibCantera.mix_updatePhases(mp.handle))
    return mp
end

# ---- sizes ------------------------------------------------------------------

"Number of phases in the mixture."
n_phases(mp::MultiPhase) = Int(check(LibCantera.mix_nPhases(mp.handle)))

"Total number of species across all phases."
n_species(mp::MultiPhase) = Int(check(LibCantera.mix_nSpecies(mp.handle)))

"Number of elements in the mixture."
n_elements(mp::MultiPhase) = Int(check(LibCantera.mix_nElements(mp.handle)))

"""
    element_index(mp::MultiPhase, name) -> Int

1-based global index of element `name`, or `0` if not present.
"""
function element_index(mp::MultiPhase, name::AbstractString)
    idx = LibCantera.mix_elementIndex(mp.handle, name)
    idx == -1 && return 0
    return Int(idx) + 1
end

# ---- scalar state -----------------------------------------------------------

"Temperature [K]."
temperature(mp::MultiPhase) = checkd(LibCantera.mix_temperature(mp.handle))

"Set the mixture temperature [K]."
function set_temperature!(mp::MultiPhase, T::Real)
    check(LibCantera.mix_setTemperature(mp.handle, Float64(T)))
    return mp
end

"Pressure [Pa]."
pressure(mp::MultiPhase) = checkd(LibCantera.mix_pressure(mp.handle))

"Set the mixture pressure [Pa]."
function set_pressure!(mp::MultiPhase, P::Real)
    check(LibCantera.mix_setPressure(mp.handle, Float64(P)))
    return mp
end

"Minimum temperature [K] of the mixture's valid range."
min_temp(mp::MultiPhase) = checkd(LibCantera.mix_minTemp(mp.handle))

"Maximum temperature [K] of the mixture's valid range."
max_temp(mp::MultiPhase) = checkd(LibCantera.mix_maxTemp(mp.handle))

"Total electric charge [C/mol] of the mixture."
charge(mp::MultiPhase) = checkd(LibCantera.mix_charge(mp.handle))

# ---- mixture thermodynamic properties ---------------------------------------

"Total enthalpy [J]."
enthalpy(mp::MultiPhase) = checkd(LibCantera.mix_enthalpy(mp.handle))

"Total entropy [J/K]."
entropy(mp::MultiPhase) = checkd(LibCantera.mix_entropy(mp.handle))

"Total Gibbs free energy [J]."
gibbs(mp::MultiPhase) = checkd(LibCantera.mix_gibbs(mp.handle))

"Total heat capacity at constant pressure [J/K]."
cp(mp::MultiPhase) = checkd(LibCantera.mix_cp(mp.handle))

"Total volume [m^3]."
volume(mp::MultiPhase) = checkd(LibCantera.mix_volume(mp.handle))

# ---- composition ------------------------------------------------------------

"Number of moles in phase `n` (1-based)."
phase_moles(mp::MultiPhase, n::Integer) =
    checkd(LibCantera.mix_phaseMoles(mp.handle, Int32(n - 1)))

"Set the number of moles in phase `n` (1-based)."
function set_phase_moles!(mp::MultiPhase, n::Integer, moles::Real)
    check(LibCantera.mix_setPhaseMoles(mp.handle, Int32(n - 1), Float64(moles)))
    return mp
end

"Moles of global species `k` (1-based)."
species_moles(mp::MultiPhase, k::Integer) =
    checkd(LibCantera.mix_speciesMoles(mp.handle, Int32(k - 1)))

"Moles of element `m` (1-based)."
element_moles(mp::MultiPhase, m::Integer) =
    checkd(LibCantera.mix_elementMoles(mp.handle, Int32(m - 1)))

"Mole fraction of global species `k` (1-based)."
mole_fraction(mp::MultiPhase, k::Integer) =
    checkd(LibCantera.mix_moleFraction(mp.handle, Int32(k - 1)))

"""
    n_atoms(mp::MultiPhase, k, m) -> Float64

Number of atoms of element `m` in global species `k` (both 1-based).
"""
n_atoms(mp::MultiPhase, k::Integer, m::Integer) =
    checkd(LibCantera.mix_nAtoms(mp.handle, Int32(k - 1), Int32(m - 1)))

"""
    set_moles!(mp::MultiPhase, moles)

Set the moles of all global species from a numeric vector, or from a
composition string (e.g. `"CH4:1, O2:2"`).
"""
function set_moles!(mp::MultiPhase, moles::AbstractVector{<:Real})
    n = as_f64(moles)
    check(LibCantera.mix_setMoles(mp.handle, Int32(length(n)), pointer(n)))
    return mp
end
function set_moles!(mp::MultiPhase, moles::AbstractString)
    check(LibCantera.mix_setMolesByName(mp.handle, moles))
    return mp
end

"""
    chemical_potentials(mp::MultiPhase) -> Vector{Float64}

Chemical potentials [J/kmol] of all global species.
"""
chemical_potentials(mp::MultiPhase) =
    get_array(n_species(mp),
        (n, b) -> LibCantera.mix_getChemPotentials(mp.handle, n, b))

"""
    equilibrate!(mp::MultiPhase, XY; solver="auto", rtol=1e-9, max_steps=1000,
                 max_iter=200, estimate_equil=0)

Bring the mixture to chemical equilibrium holding the property pair `XY`
(e.g. `"TP"`) fixed.
"""
function equilibrate!(mp::MultiPhase, XY::AbstractString; solver="auto",
                      rtol=1e-9, max_steps=1000, max_iter=200, estimate_equil=0)
    check(LibCantera.mix_equilibrate(mp.handle, XY, solver, Float64(rtol),
                                     Int32(max_steps), Int32(max_iter),
                                     Int32(estimate_equil)))
    return mp
end

function Base.show(io::IO, mp::MultiPhase)
    if mp.closed
        print(io, "MultiPhase(<closed>)")
    else
        print(io, "MultiPhase(", n_phases(mp), " phases, ", n_species(mp),
              " species)")
    end
end

# note: `cp` is available as `Cantera.cp` to avoid the clash with `Base.cp`
export MultiPhase, add_phase!, init!, update_phases!, n_phases,
       set_temperature!, set_pressure!, min_temp, max_temp, charge,
       enthalpy, entropy, gibbs, volume, phase_moles, set_phase_moles!,
       species_moles, element_moles, mole_fraction, n_atoms, set_moles!
