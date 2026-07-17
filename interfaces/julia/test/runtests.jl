# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test

# Reference values were generated with the Python Cantera interface (v3.2.0) for
# gri30.yaml and are checked in below.  See test/reference_values.md for the
# script that regenerates them.

const MECH = "gri30.yaml"

@testset "Cantera.jl" begin

    @testset "Library loading & version" begin
        v = cantera_version()
        @test v isa String
        @test !isempty(v)
        @test occursin(r"^\d+\.\d+", v)
    end

    @testset "Solution construction & finalization" begin
        gas = Solution(MECH)
        @test n_species(gas) == 53
        @test n_reactions(gas) == 325
        gas2 = Solution(MECH, "gri30")
        @test n_species(gas2) == 53
        # close! is idempotent and must not crash.
        close!(gas)
        close!(gas)
        @test occursin("closed", sprint(show, gas))
        close!(gas2)
    end

    @testset "Sub-phase views keep their Solution alive" begin
        # The views returned by `thermo`/`kinetics`/`transport` borrow handles
        # owned by the Solution, so they hold a reference to it: collecting the
        # Solution here would run its finalizer and free the handles the views
        # are still using.  Each view is built from its own Solution and is the
        # only thing referencing it, otherwise the views would keep each other's
        # owner alive and a regression in any one of them would go unnoticed.
        view(f) = let gas = Solution(MECH)
            set_TPX!(gas, 1200.0, one_atm, "CH4:1, O2:2, N2:7.52")
            f(gas)
        end

        tp = view(thermo)
        kin = view(kinetics)
        tr = view(transport)
        GC.gc(); GC.gc()

        for v in (tp, kin, tr)
            @test v.parent isa Solution
            @test !v.parent.closed
        end
        @test temperature(tp) ≈ 1200.0
        @test n_reactions(kin) == 325
        @test viscosity(tr) > 0
    end

    @testset "Error handling" begin
        @test_throws CanteraError Solution("this_file_does_not_exist.yaml")
        gas = Solution(MECH)
        # An unknown species name yields index 0 (not found), not an error.
        @test species_index(gas, "NOTASPECIES") == 0
        close!(gas)
    end

    @testset "Thermo properties vs reference" begin
        gas = Solution(MECH)
        set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")

        @test temperature(gas) ≈ 1200.0
        @test pressure(gas) ≈ one_atm rtol=1e-10
        @test density(gas) ≈ 0.2806317906177915 rtol=1e-8
        @test cp_mass(gas) ≈ 1397.2506879911464 rtol=1e-8
        @test cv_mass(gas) ≈ 1096.3671002332637 rtol=1e-8
        @test enthalpy_mass(gas) ≈ 861934.8781246373 rtol=1e-8
        @test entropy_mass(gas) ≈ 8914.227317028952 rtol=1e-8
        @test mean_molecular_weight(gas) ≈ 27.633486692015207 rtol=1e-8
        close!(gas)
    end

    @testset "Composition & species" begin
        gas = Solution(MECH)
        names = species_names(gas)
        @test length(names) == 53
        @test names[1] == "H2"
        @test "CH4" in names
        @test species_index(gas, "OH") == 5           # 1-based
        @test species_name(gas, 5) == "OH"

        mw = molecular_weights(gas)
        @test length(mw) == 53
        @test mw[species_index(gas, "CH4")] ≈ 16.043 rtol=1e-4

        # setters via string and via vector agree
        set_TPX!(gas, 300.0, one_atm, "CH4:1, O2:2")
        X = mole_fractions(gas)
        @test sum(X) ≈ 1.0 rtol=1e-10
        gas2 = Solution(MECH)
        set_TP!(gas2, 300.0, one_atm)
        set_mole_fractions!(gas2, X)
        @test mole_fractions(gas2) ≈ X rtol=1e-10

        # mass-fraction round trip
        Y = mass_fractions(gas)
        @test sum(Y) ≈ 1.0 rtol=1e-10
        close!(gas); close!(gas2)
    end

    @testset "Kinetics vs reference" begin
        gas = Solution(MECH)
        set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")

        eqs = reaction_equations(gas)
        @test length(eqs) == 325
        @test eqs[1] == "2 O + M <=> O2 + M"

        w = net_production_rates(gas)
        @test length(w) == 53
        @test w[species_index(gas, "CH4")] ≈ -8.650790742316334e-6 rtol=1e-6

        kf = forward_rate_constants(gas)
        @test length(kf) == 325
        @test kf[1] ≈ 1.0e8 rtol=1e-8

        # net = fwd - rev rates of progress
        rop_net = net_rates_of_progress(gas)
        rop_f = forward_rates_of_progress(gas)
        rop_r = reverse_rates_of_progress(gas)
        @test rop_net ≈ rop_f .- rop_r rtol=1e-10

        # net production = creation - destruction
        @test net_production_rates(gas) ≈
              creation_rates(gas) .- destruction_rates(gas) rtol=1e-8
        close!(gas)
    end

    @testset "Transport vs reference" begin
        gas = Solution(MECH)
        set_TPX!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
        @test transport_model(gas) == "mixture-averaged"
        @test viscosity(gas) ≈ 4.6872182846913704e-5 rtol=1e-6
        @test thermal_conductivity(gas) ≈ 0.08965525008638553 rtol=1e-6

        d = mix_diff_coeffs(gas)
        @test length(d) == 53
        @test all(d .> 0)
        B = binary_diff_coeffs(gas)
        @test size(B) == (53, 53)
        close!(gas)
    end

    @testset "In-place array getters" begin
        gas = Solution(MECH)
        set_TPX!(gas, 1500.0, 2 * one_atm, "CH4:1, O2:2, N2:7.52")
        out = Vector{Float64}(undef, n_species(gas))
        net_production_rates!(out, gas)
        @test out == net_production_rates(gas)
        mole_fractions!(out, gas)
        @test out ≈ mole_fractions(gas)
        close!(gas)
    end

    @testset "Reactor time integration vs reference" begin
        gas = Solution(MECH)
        set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")
        r = IdealGasReactor(gas)
        net = ReactorNet(r)
        @test time(net) ≈ 0.0
        advance!(net, 1e-3)
        @test time(net) ≈ 1e-3
        # ignition to equilibrium-ish temperature
        @test temperature(r) ≈ 2869.724616612829 rtol=1e-3
        close!(net); close!(r); close!(gas)
    end

    @testset "Analytic kinetic derivatives" begin
        gas = Solution(MECH)
        set_TPX!(gas, 1400.0, one_atm, "CH4:1, O2:2, N2:7.52")
        d = net_production_rates_ddT(gas)
        @test length(d) == 53
        @test all(isfinite, d)
        # value check against the Python Cantera reference
        @test d[species_index(gas, "CH4")] ≈ -3.3322745901816125e-6 rtol=1e-6
        # in-place variant agrees
        out = similar(d); net_production_rates_ddT!(out, gas)
        @test out == d
        close!(gas)
    end

    include("test_mechanisms.jl")
    include("test_thermo.jl")
    include("test_kinetics.jl")
    include("test_transport.jl")
    include("test_reactor.jl")
    include("test_solutionarray.jl")
    include("test_func1.jl")
    include("test_multiphase.jl")
    include("test_connectors.jl")
    include("test_rdiag.jl")
    include("test_onedim.jl")
end
