# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Tests for reactor connectors, reservoirs, surfaces and sensitivities.

using Cantera, Test

# Names not necessarily exported yet are reached through the module.
import Cantera: Reservoir, Wall, MassFlowController, Valve, PressureController,
                ReactorSurface, mass_flow_rate, set_mass_flow_rate!,
                device_coefficient, set_device_coefficient!,
                heat_rate, expansion_rate, area, set_area!,
                set_heat_transfer_coeff!, add_sensitivity_reaction!,
                n_sens_params, sensitivity, rtol, atol,
                set_sensitivity_tolerances!, connector_type, name

const CONN_MECH = "h2o2.yaml"

@testset "MassFlowController" begin
    gas = Solution(CONN_MECH)
    set_TPX!(gas, 1000.0, one_atm, "H2:2,O2:1")
    gas2 = Solution(CONN_MECH)
    set_TPX!(gas2, 1000.0, one_atm, "H2:2,O2:1")
    up = Reservoir(gas2)
    r = IdealGasReactor(gas)
    mfc = MassFlowController(up, r; mdot=0.1)
    net = ReactorNet([r])
    advance!(net, 0.0)   # initialize the network so the device is ready
    @test mass_flow_rate(mfc) ≈ 0.1
    @test device_coefficient(mfc) ≈ 0.1
    @test connector_type(mfc) == "MassFlowController"
end

@testset "Wall heat transfer" begin
    g1 = Solution(CONN_MECH); set_TPX!(g1, 1000.0, one_atm, "H2:2,O2:1")
    g2 = Solution(CONN_MECH); set_TPX!(g2,  300.0, one_atm, "H2:2,O2:1")
    r1 = IdealGasReactor(g1); r2 = IdealGasReactor(g2)
    set_initial_volume!(r1, 1.0); set_initial_volume!(r2, 1.0)
    set_chemistry_enabled!(r1, false); set_chemistry_enabled!(r2, false)
    w = Wall(r1, r2; A=1.0, U=100.0)
    @test area(w) ≈ 1.0
    net = ReactorNet([r1, r2])
    T1_0 = temperature(r1); T2_0 = temperature(r2)
    advance!(net, 1.0)
    T1 = temperature(r1); T2 = temperature(r2)
    # temperatures move toward each other
    @test T1 < T1_0
    @test T2 > T2_0
    @test T1 > T2
    q = heat_rate(w)
    @test isfinite(q) && q > 0
end

@testset "Valve" begin
    g1 = Solution(CONN_MECH); set_TPX!(g1, 500.0, 3 * one_atm, "H2:2,O2:1")
    g2 = Solution(CONN_MECH); set_TPX!(g2, 500.0, 1 * one_atm, "H2:2,O2:1")
    r1 = IdealGasReactor(g1); r2 = IdealGasReactor(g2)
    set_initial_volume!(r1, 1.0); set_initial_volume!(r2, 1.0)
    set_chemistry_enabled!(r1, false); set_chemistry_enabled!(r2, false)
    K = 1e-5
    v = Valve(r1, r2; K=K)
    @test connector_type(v) == "Valve"
    net = ReactorNet([r1, r2])
    advance!(net, 0.0)
    p1_0 = pressure(r1); p2_0 = pressure(r2)
    @test mass_flow_rate(v) ≈ K * (p1_0 - p2_0) rtol=1e-6   # flow high -> low
    @test mass_flow_rate(v) > 0
    advance!(net, 5.0)
    @test pressure(r1) < p1_0                # high side drops
    @test pressure(r2) > p2_0                # low side rises
    @test pressure(r1) ≈ pressure(r2) rtol=1e-2   # pressures equilibrate
end

@testset "PressureController" begin
    up = Solution(CONN_MECH); set_TPX!(up, 800.0, one_atm, "H2:2,O2:1")
    rg = Solution(CONN_MECH); set_TPX!(rg, 800.0, one_atm, "H2:2,O2:1")
    down = Solution(CONN_MECH); set_TPX!(down, 800.0, one_atm, "H2:2,O2:1")
    upstream = Reservoir(up)
    r = IdealGasReactor(rg); set_chemistry_enabled!(r, false)
    downstream = Reservoir(down)

    mfc = MassFlowController(upstream, r; mdot=0.05)
    pc = PressureController(r, downstream; primary=mfc, K=1e-4)
    @test connector_type(pc) == "PressureController"
    net = ReactorNet([r])
    advance!(net, 0.0)

    @test mass_flow_rate(pc) ≈ mass_flow_rate(mfc) rtol=1e-6
    @test isfinite(mass_flow_rate(pc))
end

@testset "Sensitivity" begin
    gas = Solution(CONN_MECH)
    set_TPX!(gas, 1000.0, one_atm, "H2:2,O2:1")
    r = IdealGasReactor(gas)
    net = ReactorNet([r])
    add_sensitivity_reaction!(r, 1)      # 1-based -> reaction 0
    @test n_sens_params(r) == 1
    set_sensitivity_tolerances!(net; rtol=1e-6, atol=1e-6)
    advance!(net, 2e-4)
    s = sensitivity(net, "temperature", 1, r)
    @test isfinite(s)
    @test rtol(net) > 0
    @test atol(net) > 0
end
