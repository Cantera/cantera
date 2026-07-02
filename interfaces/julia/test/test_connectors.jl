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

const MECH = "h2o2.yaml"

@testset "MassFlowController" begin
    gas = Solution(MECH)
    set_TPX!(gas, 1000.0, one_atm, "H2:2,O2:1")
    gas2 = Solution(MECH)
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
    g1 = Solution(MECH); set_TPX!(g1, 1000.0, one_atm, "H2:2,O2:1")
    g2 = Solution(MECH); set_TPX!(g2,  300.0, one_atm, "H2:2,O2:1")
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

@testset "Sensitivity" begin
    gas = Solution(MECH)
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
