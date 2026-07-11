# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test

# Reference values from Python Cantera 3.2.0, gri30 at (1200 K, 1 atm,
# CH4:1, O2:2, N2:7.52), unless noted otherwise.
@testset "ThermoPhase (Python parity)" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")

    @testset "scalar properties" begin
        @test isothermal_compressibility(gas) ≈ 9.869232667160127e-06 rtol=1e-8
        @test thermal_expansion_coeff(gas)    ≈ 0.0008333333333333334 rtol=1e-8
        @test max_temp(gas) ≈ 3000.0
        @test min_temp(gas) ≈ 300.0
        @test sound_speed(gas)        ≈ 678.3425211864376 rtol=1e-8
        @test volume_mass(gas)        ≈ 3.5633881599749238 rtol=1e-10
        @test volume_mole(gas)        ≈ 98.46883929715162 rtol=1e-10
        @test reference_pressure(gas) ≈ 101325.0
        @test electric_potential(gas) ≈ 0.0
        @test equivalence_ratio(gas)  ≈ 1.0 rtol=1e-10
    end

    @testset "elements and atoms" begin
        iC = element_index(gas, "C")
        @test elemental_mass_fraction(gas, iC) ≈ 0.04131690114779184 rtol=1e-10
        @test elemental_mole_fraction(gas, iC) ≈ 0.0415973377703827  rtol=1e-10
        @test atomic_weights(gas) ≈ [15.999, 1.008, 12.011, 14.007, 39.95] rtol=1e-8
        @test all(charges(gas) .== 0.0)
    end

    @testset "set_equivalence_ratio! round-trip" begin
        set_equivalence_ratio!(gas, 0.5, "CH4:1", "O2:1, N2:3.76")
        @test equivalence_ratio(gas) ≈ 0.5 rtol=1e-8
    end

    @testset "standard-state properties" begin
        set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
        k = species_index(gas, "CH4")
        @test standard_cp_R(gas)[k]            ≈ 9.790762829272001  rtol=1e-8
        @test standard_enthalpies_RT(gas)[k]   ≈ -2.0465605875189334 rtol=1e-8
        @test standard_entropies_R(gas)[k]     ≈ 31.561123052390712 rtol=1e-8
        @test standard_gibbs_RT(gas)[k]        ≈ -33.607683639909645 rtol=1e-8
        @test standard_int_energies_RT(gas)[k] ≈ -3.0465605875189334 rtol=1e-8
        # G/RT = H/RT - S/R identity
        @test standard_gibbs_RT(gas)[k] ≈
              standard_enthalpies_RT(gas)[k] - standard_entropies_R(gas)[k] rtol=1e-10
    end

    close!(gas)
end


@testset "two-property state setters (round-trip)" begin
    gas = Solution("gri30.yaml")
    comp = "CH4:1.0, O2:2.0, N2:7.52"
    set_TPX!(gas, 1200.0, one_atm, comp)
    T0 = temperature(gas); p0 = pressure(gas); rho0 = density(gas)
    h0 = enthalpy_mass(gas); u0 = internal_energy_mass(gas)
    s0 = entropy_mass(gas); v0 = volume_mass(gas)

    perturb() = set_TPX!(gas, 600.0, 2 * one_atm, comp)  # composition unchanged

    @testset "set_HP!" begin
        perturb(); set_HP!(gas, h0, p0)
        @test temperature(gas) ≈ T0 rtol=1e-6
        @test pressure(gas)    ≈ p0 rtol=1e-10
    end
    @testset "set_UV!" begin
        perturb(); set_UV!(gas, u0, v0)
        @test temperature(gas) ≈ T0   rtol=1e-6
        @test density(gas)     ≈ rho0 rtol=1e-8
    end
    @testset "set_SP!" begin
        perturb(); set_SP!(gas, s0, p0)
        @test temperature(gas) ≈ T0 rtol=1e-6
        @test pressure(gas)    ≈ p0 rtol=1e-10
    end
    @testset "set_SV!" begin
        perturb(); set_SV!(gas, s0, v0)
        @test temperature(gas) ≈ T0   rtol=1e-6
        @test density(gas)     ≈ rho0 rtol=1e-8
    end
    @testset "set_DP!" begin
        perturb(); set_DP!(gas, rho0, p0)
        @test density(gas)  ≈ rho0 rtol=1e-10
        @test pressure(gas) ≈ p0   rtol=1e-10
    end
    @testset "set_TD!" begin
        perturb(); set_TD!(gas, T0, rho0)
        @test temperature(gas) ≈ T0   rtol=1e-10
        @test density(gas)     ≈ rho0 rtol=1e-10
    end

    close!(gas)
end

@testset "equilibrate! modes" begin
    @testset "TP is a fixed point" begin
        gas = Solution("gri30.yaml")
        set_TPX!(gas, 2000.0, one_atm, "CH4:1, O2:2, N2:7.52")
        equilibrate!(gas, "TP")
        @test temperature(gas) ≈ 2000.0 rtol=1e-8
        @test pressure(gas)    ≈ one_atm rtol=1e-8
        X1 = mole_fractions(gas)
        equilibrate!(gas, "TP")          # idempotent
        @test mole_fractions(gas) ≈ X1 atol=1e-8
        close!(gas)
    end
    @testset "HP conserves enthalpy & pressure (adiabatic combustion)" begin
        gas = Solution("gri30.yaml")
        set_TPX!(gas, 1500.0, one_atm, "H2:2, O2:1, N2:4")
        h0 = enthalpy_mass(gas)
        equilibrate!(gas, "HP")
        @test enthalpy_mass(gas) ≈ h0      rtol=1e-6
        @test pressure(gas)      ≈ one_atm rtol=1e-8
        @test temperature(gas) > 1500.0    # exothermic: burned gas is hotter
        close!(gas)
    end
    @testset "UV conserves internal energy & volume" begin
        gas = Solution("gri30.yaml")
        set_TPX!(gas, 1500.0, one_atm, "H2:2, O2:1, N2:4")
        u0 = internal_energy_mass(gas); v0 = volume_mass(gas)
        equilibrate!(gas, "UV")
        @test internal_energy_mass(gas) ≈ u0 rtol=1e-6
        @test volume_mass(gas)          ≈ v0 rtol=1e-8
        close!(gas)
    end
end
