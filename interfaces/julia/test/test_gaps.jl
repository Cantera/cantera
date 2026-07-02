using Cantera
using Test
using LinearAlgebra
import Libdl

# Features added to close remaining Python/MATLAB parity gaps: the composition
# Jacobian (_ddX, via a densifying CLib recipe) and ReactorNet state/component
# introspection.  Values checked against Python Cantera 3.x for gri30.yaml.
#
# The _ddX and ReactorNet-state getters are CLib functions added on this branch;
# they are only present in a libcantera built from this Cantera. Against an older
# library the symbol is absent, so those checks are skipped gracefully.
_have(sym) = Libdl.dlsym_e(Libdl.dlopen(Cantera.LibCantera.libcantera[]), sym) != C_NULL

if !_have(:kin_getNetProductionRates_ddX)
    @info "Skipping _ddX tests: linked libcantera lacks the ctkin _ddX recipes " *
          "(build Cantera from this branch to enable)."
    @testset "Composition Jacobian (_ddX)" begin
        @test_skip false
    end
else
@testset "Composition Jacobian (_ddX)" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1400.0, one_atm, "CH4:1, O2:2, N2:7.52")

    J = net_production_rates_ddX(gas)
    @test size(J) == (53, 53)
    @test J[species_index(gas, "CH4"), species_index(gas, "CH4")] ≈
          -0.0022867583538256886 rtol=1e-8
    @test norm(J) ≈ 7289124.658794 rtol=1e-8

    R = net_rates_of_progress_ddX(gas)
    @test size(R) == (n_reactions(gas), n_species(gas))
    @test all(isfinite, R)
    close!(gas)
end
end  # _have(:kin_getNetProductionRates_ddX)

if !_have(:reactornet_neq)
    @info "Skipping ReactorNet state tests: linked libcantera lacks the " *
          "ctreactornet state recipes (build Cantera from this branch to enable)."
    @testset "ReactorNet state & components" begin
        @test_skip false
    end
else
@testset "ReactorNet state & components" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")
    r = IdealGasReactor(gas)
    net = ReactorNet(r)

    @test n_components(net) == 3 + n_species(gas)      # mass, volume, temperature + species
    y = state(net)
    @test length(y) == n_components(net)
    @test all(isfinite, y)

    names = component_names(net)
    @test length(names) == n_components(net)
    @test occursin("mass", names[1])
    @test occursin("temperature", names[3])
    close!(net); close!(r); close!(gas)
end
end  # _have(:reactornet_neq)

@testset "SolutionArray" begin
    gas = Solution("gri30.yaml")
    # snapshot a short reactor history
    set_TPX!(gas, 1200.0, one_atm, "H2:2, O2:1, N2:4")
    r = IdealGasReactor(gas); net = ReactorNet(r)
    states = SolutionArray(gas)
    for t in range(0, 1e-3; length=5)
        advance!(net, t)
        append!(states, gas)
    end
    @test length(states) == 5
    T = temperature(states)
    @test length(T) == 5
    @test T[end] > T[1] + 500               # constant-volume ignition heats up
    @test maximum(T) > 2000
    @test size(mass_fractions(states)) == (n_species(gas), 5)
    @test all(density(states) .> 0)
    @test extract(states, enthalpy_mass) isa Vector{Float64}

    # fixed-size construction + set_state! + CSV round-trip
    sa = SolutionArray(gas, 3)
    @test length(sa) == 3
    set_state!(sa, 1; T=300.0, P=one_atm, X="CH4:1, O2:2")
    @test temperature(sa)[1] ≈ 300.0
    path = tempname() * ".csv"
    write_csv(sa, path)
    @test isfile(path)
    @test occursin("T,P,", readline(path))
    rm(path; force=true)
    close!(net); close!(r); close!(gas)
end
