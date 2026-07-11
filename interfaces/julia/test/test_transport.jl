# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test
using LinearAlgebra: diag


@testset "Transport / mixture-averaged" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1200.0, one_atm, "CH4:1, O2:2, N2:7.52")
    nsp = n_species(gas)

    @test viscosity(gas) > 0
    @test thermal_conductivity(gas) > 0

    D = mix_diff_coeffs(gas)
    @test length(D) == nsp
    @test all(D .> 0)

    Dbin = binary_diff_coeffs(gas)
    @test size(Dbin) == (nsp, nsp)
    @test all(isfinite, Dbin)
    @test Dbin ≈ transpose(Dbin)          # binary coefficients are symmetric
    @test all(diag(Dbin) .> 0)

    # Diffusion coefficients grow with temperature at fixed pressure.
    D1200 = mix_diff_coeffs(gas)
    set_TPX!(gas, 1800.0, one_atm, "CH4:1, O2:2, N2:7.52")
    @test all(mix_diff_coeffs(gas) .> D1200)

    close!(gas)
end

@testset "Transport / multicomponent" begin
    gas = Solution("gri30.yaml"; transport="multicomponent")
    set_TPX!(gas, 1200.0, one_atm, "CH4:1, O2:2, N2:7.52")
    nsp = n_species(gas)

    @test transport_model(gas) == "multicomponent"
    @test viscosity(gas) > 0
    @test thermal_conductivity(gas) > 0

    Dmulti = multi_diff_coeffs(gas)
    @test size(Dmulti) == (nsp, nsp)
    @test all(isfinite, Dmulti)           # multicomponent coeffs may be negative

    DT = thermal_diff_coeffs(gas)
    @test length(DT) == nsp
    @test all(isfinite, DT)
    @test any(!=(0.0), DT)

    close!(gas)
end
