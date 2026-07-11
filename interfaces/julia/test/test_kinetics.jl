# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test
using LinearAlgebra

# Reference values from Python Cantera 3.2.0, gri30 at (1200 K, 1 atm,
# CH4:1, O2:2, N2:7.52), unless noted otherwise.
@testset "Kinetics (Python parity)" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")

    @testset "stoichiometry and rates" begin
        kCH4 = species_index(gas, "CH4")
        @test sum(reactant_stoich_coeffs(gas)[kCH4, :]) ≈ 6.0
        @test size(product_stoich_coeffs(gas)) == (n_species(gas), n_reactions(gas))
        @test delta_standard_enthalpy(gas)[1] ≈ -506662565.3788193 rtol=1e-8
        @test forward_rates_of_progress_ddT(gas)[1] ≈ 0.0 atol=1e-30
        @test multiplier(gas, 1) ≈ 1.0
        set_multiplier!(gas, 1, 2.0); @test multiplier(gas, 1) ≈ 2.0
    end

    @testset "heat release / production rates" begin
        @test heat_release_rate(gas) ≈ -2068.555280246291 rtol=1e-8
        hpr = heat_production_rates(gas)
        @test length(hpr) == n_reactions(gas)
        # sum over reactions equals the total heat release rate
        @test sum(hpr) ≈ heat_release_rate(gas) rtol=1e-10
        @test sum(hpr) ≈ -2068.5552802462907 rtol=1e-8
    end

    close!(gas)
end

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
