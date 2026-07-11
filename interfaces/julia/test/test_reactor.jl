# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test

@testset "ReactorNet state & components" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")
    r = IdealGasReactor(gas)
    net = ReactorNet(r)

    # mass, volume, temperature + species
    @test n_components(net) == 3 + n_species(gas)
    y = state(net)
    @test length(y) == n_components(net)
    @test all(isfinite, y)

    names = component_names(net)
    @test length(names) == n_components(net)
    @test occursin("mass", names[1])
    @test occursin("temperature", names[3])
    close!(net); close!(r); close!(gas)
end

@testset "ConstPressureReactor ignition" begin
    for ctor in (ConstPressureReactor, IdealGasConstPressureReactor)
        gas = Solution("gri30.yaml")
        set_TPX!(gas, 1200.0, one_atm, "H2:2, O2:1, N2:4")
        p0 = pressure(gas)

        # Independent target: HP equilibrium of the same initial mixture.
        ad = Solution("gri30.yaml")
        set_TPX!(ad, 1200.0, one_atm, "H2:2, O2:1, N2:4")
        equilibrate!(ad, "HP")
        Tad = temperature(ad)

        r = ctor(gas)
        net = ReactorNet([r])
        advance!(net, 1.0)               # well past ignition

        @test pressure(r)    ≈ p0  rtol=1e-4    # constant pressure
        @test temperature(r) ≈ Tad rtol=2e-2    # adiabatic flame temperature
        close!(net); close!(r); close!(ad); close!(gas)
    end
end

@testset "Reservoir holds state" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 800.0, one_atm, "CH4:1, O2:2, N2:7.52")
    res = Reservoir(gas)
    T0 = temperature(res); p0 = pressure(res)

    driver = Solution("gri30.yaml")
    set_TPX!(driver, 1200.0, one_atm, "H2:2, O2:1, N2:4")
    r = IdealGasReactor(driver)
    MassFlowController(res, r; mdot=0.01)   # reservoir feeds the reactor
    net = ReactorNet([r])
    advance!(net, 0.5)

    @test temperature(res) ≈ T0 rtol=1e-10   # unchanged
    @test pressure(res)    ≈ p0 rtol=1e-10
    close!(net); close!(r); close!(res); close!(gas); close!(driver)
end
