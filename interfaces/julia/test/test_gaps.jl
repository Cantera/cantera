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

# Scalar thermo/kinetics getters that fill the last parity gaps vs. Python.
# Reference values from Python Cantera 3.2.0, gri30 at (1200 K, 1 atm,
# CH4:1, O2:2, N2:7.52).
@testset "Extra scalar getters (Python parity)" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
    @test isothermal_compressibility(gas) ≈ 9.869232667160127e-06 rtol=1e-8
    @test thermal_expansion_coeff(gas)    ≈ 0.0008333333333333334 rtol=1e-8
    @test max_temp(gas) ≈ 3000.0
    @test min_temp(gas) ≈ 300.0
    @test heat_release_rate(gas) ≈ -2068.555280246291 rtol=1e-8
    close!(gas)
end

# Wrapped/derived getters closing the remaining Python parity gaps.
# Reference values from Python Cantera 3.2.0, gri30 at (1200 K, 1 atm,
# CH4:1, O2:2, N2:7.52).
@testset "Wrapped parity getters" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
    @test sound_speed(gas)        ≈ 678.3425211864376 rtol=1e-8
    @test volume_mass(gas)        ≈ 3.5633881599749238 rtol=1e-10
    @test volume_mole(gas)        ≈ 98.46883929715162 rtol=1e-10
    @test reference_pressure(gas) ≈ 101325.0
    @test electric_potential(gas) ≈ 0.0
    @test equivalence_ratio(gas)  ≈ 1.0 rtol=1e-10
    iC = element_index(gas, "C")
    @test elemental_mass_fraction(gas, iC) ≈ 0.04131690114779184 rtol=1e-10
    @test elemental_mole_fraction(gas, iC) ≈ 0.0415973377703827  rtol=1e-10
    @test atomic_weights(gas) ≈ [15.999, 1.008, 12.011, 14.007, 39.95] rtol=1e-8
    @test all(charges(gas) .== 0.0)
    # kinetics
    kCH4 = species_index(gas, "CH4")
    @test sum(reactant_stoich_coeffs(gas)[kCH4, :]) ≈ 6.0
    @test size(product_stoich_coeffs(gas)) == (n_species(gas), n_reactions(gas))
    @test delta_standard_enthalpy(gas)[1] ≈ -506662565.3788193 rtol=1e-8
    @test forward_rates_of_progress_ddT(gas)[1] ≈ 0.0 atol=1e-30
    @test multiplier(gas, 1) ≈ 1.0
    set_multiplier!(gas, 1, 2.0); @test multiplier(gas, 1) ≈ 2.0
    # set_equivalence_ratio! round-trip
    set_equivalence_ratio!(gas, 0.5, "CH4:1", "O2:1, N2:3.76")
    @test equivalence_ratio(gas) ≈ 0.5 rtol=1e-8
    close!(gas)
end

# Standard-state properties are CLib functions added on this branch (ctthermo
# recipes); skip against a libcantera that predates them.
if !_have(:thermo_getCp_R)
    @info "Skipping standard-state tests: linked libcantera lacks ctthermo standard-state recipes."
else
@testset "Standard-state properties (Python parity)" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
    k = species_index(gas, "CH4")
    @test standard_cp_R(gas)[k]            ≈ 9.790762829272001  rtol=1e-8
    @test standard_enthalpies_RT(gas)[k]   ≈ -2.0465605875189334 rtol=1e-8
    @test standard_entropies_R(gas)[k]     ≈ 31.561123052390712 rtol=1e-8
    @test standard_gibbs_RT(gas)[k]        ≈ -33.607683639909645 rtol=1e-8
    @test standard_int_energies_RT(gas)[k] ≈ -3.0465605875189334 rtol=1e-8
    # G/RT = H/RT - S/R identity
    @test standard_gibbs_RT(gas)[k] ≈ standard_enthalpies_RT(gas)[k] - standard_entropies_R(gas)[k] rtol=1e-10
    close!(gas)
end
end  # _have(:thermo_getCp_R)

@testset "heat_production_rates" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
    hpr = heat_production_rates(gas)
    @test length(hpr) == n_reactions(gas)
    # sum over reactions equals the total heat release rate
    @test sum(hpr) ≈ heat_release_rate(gas) rtol=1e-10
    @test sum(hpr) ≈ -2068.5552802462907 rtol=1e-8
    close!(gas)
end
