using Cantera
using Test
import Libdl

# Reference computed with Python `cantera.FreeFlame` (gri30, stoichiometric
# CH4/air at 300 K, 1 atm, width=0.03, refine ratio=3/slope=0.06/curve=0.12,
# solve auto=True):
#   python3 -c "import cantera as ct; g=ct.Solution('gri30.yaml');
#     g.TPX=300,ct.one_atm,'CH4:1,O2:2,N2:7.52';
#     f=ct.FreeFlame(g,width=0.03); f.set_refine_criteria(ratio=3,slope=0.06,curve=0.12);
#     f.solve(loglevel=0,auto=True); print(repr(f.velocity[0]))"
#   -> Su = 0.3809265424901588 m/s, final grid = 196 points spanning [0, 0.06] m,
#      Tmax = 2230.74 K.
#
# The full `auto=true` solve reproduces this to ~6 significant figures but takes
# ~200 s (it faithfully mirrors Python's staged multi-grid + domain-widening
# schedule to land on the identical 196-point grid). To keep the default test
# suite fast, that exact-match check runs only when CANTERA_EXACT_FLAME is set;
# otherwise a quick coarse-grid smoke test exercises the same code paths.
const CT = Cantera
const SU_REF = 0.3809265424901588
const NPTS_REF = 196

# The 1-D flame façade uses `ctdomain`/`ctonedim` functions as they exist on this
# branch (e.g. `domain_domainType`). An older libcantera may lack or have renamed
# them; in that case skip the flame tests gracefully.
_have_onedim() =
    Libdl.dlsym_e(Libdl.dlopen(CT.LibCantera.libcantera[]), :domain_domainType) != C_NULL

if !_have_onedim()
    @info "Skipping 1-D flame tests: linked libcantera lacks the current " *
          "ctdomain/ctonedim API (build Cantera from this branch to enable)."
    @testset "OneDim / FreeFlame" begin
        @test_skip false
    end
else
@testset "OneDim / FreeFlame" begin
    gas = CT.Solution("gri30.yaml")
    CT.set_TPX!(gas, 300.0, 101325.0, "CH4:1, O2:2, N2:7.52")

    flame = CT.FreeFlame(gas; width=0.03)
    npts0 = CT.n_points(flame)

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

    if haskey(ENV, "CANTERA_EXACT_FLAME")
        # Full staged auto-solve: reproduces Python's grid and flame speed.
        t0 = time()
        CT.solve!(flame; loglevel=0, auto=true)
        dt = time() - t0
        Su = CT.flame_speed(flame)
        z = CT.grid(flame)
        Tmax = maximum(CT.flame_T(flame))
        err = abs(Su - SU_REF) / SU_REF
        println("Su=$(Su) m/s  ref=$(SU_REF)  rel_err=$(round(err*100; digits=4))%  " *
                "Tmax=$(round(Tmax; digits=2)) K  npts=$(length(z)) (ref $(NPTS_REF))  " *
                "z_end=$(round(z[end]; digits=4)) m  solve=$(round(dt; digits=1)) s")

        @testset "flame speed matches Python (exact)" begin
            @test isapprox(Su, SU_REF; rtol=1e-3)   # agrees to ~6 sig figs in practice
            @test length(z) == NPTS_REF             # identical final grid
            @test z[end] ≈ 0.06 atol = 1e-4
            @test Tmax > 1800.0
            @test all(diff(z) .> 0)
        end
    else
        # Fast smoke test: a bounded coarse solve on the initial grid. Exercises
        # the guess/energy/solve paths and asserts a physical flame structure
        # without the expensive grid refinement. Set CANTERA_EXACT_FLAME=1 to run
        # the full solve and check the flame speed against Python.
        CT.solve!(flame; loglevel=0, auto=false, refine_grid=false)
        z = CT.grid(flame)
        T = CT.flame_T(flame)

        @testset "coarse solve (smoke)" begin
            @test all(diff(z) .> 0)                 # strictly increasing grid
            @test z[1] ≈ 0.0 atol = 1e-9
            @test length(T) == length(z)
            @test T[1] < 400.0                      # cold reactants
            @test maximum(T) > 1800.0               # hot products
            Xch4 = CT.flame_X(flame, "CH4")
            Xco2 = CT.flame_X(flame, "CO2")
            @test Xch4[1] > Xch4[end]               # fuel consumed
            @test Xco2[end] > Xco2[1]               # products formed
            @test CT.flame_speed(flame) > 0.0
        end
    end

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

    # Bounded coarse solve only (no auto/refinement): a few seconds.
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
end  # _have_onedim()
