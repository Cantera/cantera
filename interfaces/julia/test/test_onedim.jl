# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test


const CT = Cantera

@testset "OneDim / FreeFlame" begin
    gas = CT.Solution("gri30.yaml")
    CT.set_TPX!(gas, 300.0, 101325.0, "CH4:1, O2:2, N2:7.52")

    flame = CT.FreeFlame(gas; width=0.03)
    npts0 = CT.n_points(flame)
    rho_u = flame.rho_u   # unburned density [kg/m^3]

    @testset "construction" begin
        @test CT.domain_type(flame.inlet) == "inlet"
        @test CT.domain_type(flame.flow) == "free-flow"
        @test CT.domain_type(flame.outlet) == "outlet"
        comps = CT.component_names(flame.flow)
        @test "T" in comps
        @test "velocity" in comps
        # Python's initial FreeFlame grid is 8 points; setting the fixed
        # temperature inserts one anchor node, giving 9.
        @test npts0 in (8, 9)
        @test occursin("FreeFlame", sprint(show, flame))
    end

    CT.set_refine_criteria!(flame; ratio=3.0, slope=0.06, curve=0.12, prune=0.0)

    t0 = time()
    CT.solve!(flame; loglevel=0, auto=true)
    dt = time() - t0
    Su = CT.flame_speed(flame)
    z = CT.grid(flame)
    T = CT.flame_T(flame)
    Tmax = maximum(T)

    ad = CT.Solution("gri30.yaml")
    CT.set_TPX!(ad, 300.0, 101325.0, "CH4:1, O2:2, N2:7.52")
    CT.equilibrate!(ad, "HP")
    Tad = CT.temperature(ad)
    rho_b = CT.density(ad)   # burned density at the adiabatic state

    println("Su=$(round(Su; digits=4)) m/s  Tmax=$(round(Tmax; digits=1)) K  " *
            "Tad=$(round(Tad; digits=1)) K  npts=$(length(z))  " *
            "z_end=$(round(z[end]; digits=4)) m  solve=$(round(dt; digits=1)) s")

    @testset "flame is physical (behavioral)" begin
        # Burned gas reaches the adiabatic flame temperature (Python: rel=2e-2).
        @test isapprox(Tmax, Tad; rtol=2e-2)
        # Laminar flame speed in the physically expected range (Python asserts
        # flame speed to rel=1e-1 for its mixtures).
        @test isapprox(Su, 0.38; rtol=0.15)
        # Mass flux ρu is conserved across the flame: ρ_b·u_b ≈ ρ_u·S_u.
        u_b = CT.flame_velocity(flame)[end]
        @test isapprox(rho_b * u_b, rho_u * Su; rtol=0.05)
        # Grid is a valid, refined, strictly increasing mesh.
        @test all(diff(z) .> 0)
        @test length(z) > npts0
        # Cold reactants, hot products.
        @test T[1] < 400.0
        @test Tmax > 1800.0
        # Fuel consumed, products formed across the front.
        Xch4 = CT.flame_X(flame, "CH4")
        Xco2 = CT.flame_X(flame, "CO2")
        @test Xch4[1] > Xch4[end]
        @test Xco2[end] > Xco2[1]
    end

    CT.close!(ad)
    CT.close!(flame)
    @test occursin("closed", sprint(show, flame))
    CT.close!(gas)
end

@testset "OneDim / BurnerFlame" begin
    gas = CT.Solution("gri30.yaml")
    CT.set_TPX!(gas, 300.0, 101325.0, "CH4:0.9, O2:2, N2:7.52")  # lean CH4/air

    mdot = 0.06  # kg/m^2/s (burner-stabilized flat flame)
    bf = CT.BurnerFlame(gas; width=0.03, mdot=mdot)

    @testset "construction" begin
        @test CT.domain_type(bf.burner) == "inlet"
        @test CT.domain_type(bf.flow) == "unstrained-flow"
        @test CT.domain_type(bf.outlet) == "outlet"
        comps = CT.component_names(bf.flow)
        @test "T" in comps
        @test "velocity" in comps
        @test CT.n_points(bf) == 7          # [0,0.1,0.2,0.3,0.5,0.7,1.0]*width
        @test CT.burner_mdot(bf) ≈ mdot
        @test occursin("BurnerFlame", sprint(show, bf))
    end

    CT.solve!(bf; loglevel=0, auto=false, refine_grid=false)
    z = CT.grid(bf)
    T = CT.flame_T(bf)

    @testset "coarse solve (smoke)" begin
        @test all(diff(z) .> 0)                 # strictly increasing grid
        @test z[1] ≈ 0.0 atol = 1e-9
        @test length(T) == length(z)
        @test all(isfinite, T)
        @test T[1] < 400.0                      # cold burner side
        @test maximum(T) > 1500.0               # hot downstream
        Xch4 = CT.flame_X(bf, "CH4")
        Xco2 = CT.flame_X(bf, "CO2")
        @test all(isfinite, Xch4)
        @test all(isfinite, Xco2)
        @test Xch4[1] > Xch4[end]               # fuel consumed downstream
        @test Xco2[end] > Xco2[1]               # products formed downstream
    end

    CT.close!(bf)
    @test occursin("closed", sprint(show, bf))
    CT.close!(gas)
end
