# SolutionArray: a lightweight container for a batch of thermodynamic states of
# a single phase, mirroring the common uses of Python's `cantera.SolutionArray`
# (storing a reactor-integration history, a parameter sweep, or a flame profile
# and reading properties or writing CSV).
#
# Unlike the Python version this is a pure-Julia convenience: states are stored
# as (T, P, mass-fractions) and properties are evaluated by restoring each state
# into the backing `Solution` on demand.  It needs no CLib support.

"""
    SolutionArray(gas::Solution, n::Integer)
    SolutionArray(gas::Solution)

A collection of `n` thermodynamic states of `gas` (or an empty, growable array).
States are stored as temperature, pressure, and mass fractions, and are restored
into `gas` on demand to evaluate properties.

```julia
states = SolutionArray(gas)
for t in times
    advance!(net, t)
    append!(states, gas)          # snapshot the current gas state
end
temperature(states)               # Vector of T for every stored state
write_csv(states, "history.csv")
```
"""
mutable struct SolutionArray
    gas::Solution
    T::Vector{Float64}
    P::Vector{Float64}
    Y::Matrix{Float64}            # nSpecies x nStates
end

function SolutionArray(gas::Solution, n::Integer)
    nsp = n_species(gas)
    T = fill(temperature(gas), n)
    P = fill(pressure(gas), n)
    Y = repeat(mass_fractions(gas), 1, n)
    return SolutionArray(gas, T, P, Y)
end

SolutionArray(gas::Solution) =
    SolutionArray(gas, Float64[], Float64[], Matrix{Float64}(undef, n_species(gas), 0))

Base.length(sa::SolutionArray) = length(sa.T)
Base.size(sa::SolutionArray) = (length(sa),)
Base.isempty(sa::SolutionArray) = length(sa) == 0

"Restore stored state `i` (1-based) into the backing `Solution` and return it."
function restore!(sa::SolutionArray, i::Integer)
    set_TPY!(sa.gas, sa.T[i], sa.P[i], @view sa.Y[:, i])
    return sa.gas
end

"""
    append!(sa::SolutionArray, gas::Solution)

Append a snapshot of the current state of `gas` to `sa`.
"""
function Base.append!(sa::SolutionArray, gas::Solution)
    push!(sa.T, temperature(gas))
    push!(sa.P, pressure(gas))
    sa.Y = hcat(sa.Y, mass_fractions(gas))
    return sa
end

"""
    set_state!(sa, i; T, P, X=nothing, Y=nothing)

Set stored state `i` (1-based). Provide composition as mole fractions `X` or
mass fractions `Y` (string or vector); if neither is given the composition is
left unchanged.
"""
function set_state!(sa::SolutionArray, i::Integer; T, P, X=nothing, Y=nothing)
    if X !== nothing
        set_TPX!(sa.gas, T, P, X)
    elseif Y !== nothing
        set_TPY!(sa.gas, T, P, Y)
    else
        set_TP!(sa.gas, T, P)
    end
    sa.T[i] = temperature(sa.gas)
    sa.P[i] = pressure(sa.gas)
    sa.Y[:, i] = mass_fractions(sa.gas)
    return sa
end

# ---- vectorized property access --------------------------------------------

"Temperatures [K] of all stored states."
temperature(sa::SolutionArray) = copy(sa.T)

"Pressures [Pa] of all stored states."
pressure(sa::SolutionArray) = copy(sa.P)

"Mass fractions of all stored states, `nSpecies x nStates`."
mass_fractions(sa::SolutionArray) = copy(sa.Y)

"""
    extract(sa, f) -> Vector

Evaluate the scalar property function `f` (e.g. `density`, `enthalpy_mass`) for
every stored state and return the results as a vector.
"""
function extract(sa::SolutionArray, f)
    return [f(restore!(sa, i)) for i in 1:length(sa)]
end

"Mass densities [kg/m^3] of all stored states."
density(sa::SolutionArray) = extract(sa, density)

"Mole fractions of all stored states, `nSpecies x nStates`."
function mole_fractions(sa::SolutionArray)
    nsp = n_species(sa.gas)
    X = Matrix{Float64}(undef, nsp, length(sa))
    for i in 1:length(sa)
        X[:, i] = mole_fractions(restore!(sa, i))
    end
    return X
end

"""
    write_csv(sa::SolutionArray, path)

Write the stored states to a CSV file with columns `T`, `P`, and one column per
species mass fraction.
"""
function write_csv(sa::SolutionArray, path::AbstractString)
    names = species_names(sa.gas)
    open(path, "w") do io
        println(io, "T,P,", join(names, ","))
        for i in 1:length(sa)
            print(io, sa.T[i], ",", sa.P[i])
            for k in 1:size(sa.Y, 1)
                print(io, ",", sa.Y[k, i])
            end
            println(io)
        end
    end
    return path
end

Base.show(io::IO, sa::SolutionArray) =
    print(io, "SolutionArray(", length(sa), " states, ", n_species(sa.gas), " species)")
