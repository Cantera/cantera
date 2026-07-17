# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test

@testset "SolutionArray" begin
    gas = Solution("gri30.yaml")
    # snapshot a short reactor history
    set_TPX!(gas, 1200.0, one_atm, "H2:2, O2:1, N2:4")
    r = IdealGasReactor(gas); net = ReactorNet(r)
    states = SolutionArray(gas)
    for t in range(0, 1e-3; length=5)
        advance!(net, t)
        append!(states, reactor_phase(r))
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
