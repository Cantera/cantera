# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

using Cantera
using Test

# Allocations of one call, measured inside a function after a warm-up call so that
# neither compilation nor global-scope boxing is counted.
function allocations(f!, out, obj)
    f!(out, obj)
    return @allocated f!(out, obj)
end

# The contract shared by every in-place getter `f!(out, obj)`, checked against its
# allocating counterpart `f(obj)`.
function check_inplace(@nospecialize(f), @nospecialize(f!), obj)
    ref = f(obj)
    out = similar(ref)
    @test f!(out, obj) === out
    @test out == ref
    @test allocations(f!, out, obj) == 0
    @inferred f!(out, obj)

    # Wrong sizes are rejected rather than partially filled.
    if ref isa AbstractMatrix
        nr, nc = size(ref)
        @test_throws DimensionMismatch f!(similar(ref, nr + 1, nc), obj)
        @test_throws DimensionMismatch f!(similar(ref, nr, nc - 1), obj)
    else
        @test_throws DimensionMismatch f!(similar(ref, length(ref) + 1), obj)
        @test_throws DimensionMismatch f!(similar(ref, length(ref) - 1), obj)
    end

    # A contiguous view, here a column of a larger matrix, is filled in place.
    M = zeros(length(ref), 3)
    col = ref isa AbstractMatrix ? reshape(view(M, :, 2), size(ref)) : view(M, :, 2)
    f!(col, obj)
    @test vec(M[:, 2]) == vec(ref)
    @test all(iszero, M[:, 1]) && all(iszero, M[:, 3])

    # A strided view cannot be handed to the CLib.
    if ref isa AbstractVector
        @test_throws ArgumentError f!(view(zeros(2 * length(ref)), 1:2:2 * length(ref)), obj)
    end
    return nothing
end

@testset "In-place getters" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1500.0, one_atm, "CH4:1, O2:2, N2:7.52")

    @testset "thermo: $f" for (f, f!) in (
            (molecular_weights, molecular_weights!),
            (mole_fractions, mole_fractions!),
            (mass_fractions, mass_fractions!),
            (concentrations, concentrations!),
            (standard_enthalpies_RT, standard_enthalpies_RT!),
            (standard_entropies_R, standard_entropies_R!),
            (standard_gibbs_RT, standard_gibbs_RT!),
            (standard_int_energies_RT, standard_int_energies_RT!),
            (standard_cp_R, standard_cp_R!),
            (partial_molar_enthalpies, partial_molar_enthalpies!),
            (partial_molar_entropies, partial_molar_entropies!),
            (partial_molar_int_energies, partial_molar_int_energies!),
            (partial_molar_cp, partial_molar_cp!),
            (partial_molar_volumes, partial_molar_volumes!),
            (chemical_potentials, chemical_potentials!),
            (electrochemical_potentials, electrochemical_potentials!),
        )
        check_inplace(f, f!, gas)
    end

    @testset "kinetics: $f" for (f, f!) in (
            (forward_rates_of_progress, forward_rates_of_progress!),
            (reverse_rates_of_progress, reverse_rates_of_progress!),
            (net_rates_of_progress, net_rates_of_progress!),
            (equilibrium_constants, equilibrium_constants!),
            (forward_rate_constants, forward_rate_constants!),
            (reverse_rate_constants, reverse_rate_constants!),
            (delta_enthalpy, delta_enthalpy!),
            (delta_gibbs, delta_gibbs!),
            (delta_entropy, delta_entropy!),
            (delta_standard_enthalpy, delta_standard_enthalpy!),
            (delta_standard_gibbs, delta_standard_gibbs!),
            (delta_standard_entropy, delta_standard_entropy!),
            (net_production_rates, net_production_rates!),
            (creation_rates, creation_rates!),
            (destruction_rates, destruction_rates!),
            (net_production_rates_ddT, net_production_rates_ddT!),
            (net_production_rates_ddP, net_production_rates_ddP!),
            (net_production_rates_ddC, net_production_rates_ddC!),
            (creation_rates_ddT, creation_rates_ddT!),
            (creation_rates_ddP, creation_rates_ddP!),
            (creation_rates_ddC, creation_rates_ddC!),
            (destruction_rates_ddT, destruction_rates_ddT!),
            (destruction_rates_ddP, destruction_rates_ddP!),
            (destruction_rates_ddC, destruction_rates_ddC!),
            (net_rates_of_progress_ddT, net_rates_of_progress_ddT!),
            (net_rates_of_progress_ddP, net_rates_of_progress_ddP!),
            (net_rates_of_progress_ddC, net_rates_of_progress_ddC!),
            (forward_rates_of_progress_ddT, forward_rates_of_progress_ddT!),
            (forward_rates_of_progress_ddP, forward_rates_of_progress_ddP!),
            (forward_rates_of_progress_ddC, forward_rates_of_progress_ddC!),
            (reverse_rates_of_progress_ddT, reverse_rates_of_progress_ddT!),
            (reverse_rates_of_progress_ddP, reverse_rates_of_progress_ddP!),
            (reverse_rates_of_progress_ddC, reverse_rates_of_progress_ddC!),
            (net_production_rates_ddX, net_production_rates_ddX!),
            (net_rates_of_progress_ddX, net_rates_of_progress_ddX!),
        )
        check_inplace(f, f!, gas)
    end

    @testset "kinetics: keyword arguments" begin
        ref = reverse_rate_constants(gas; do_irreversible=true)
        out = similar(ref)
        @test reverse_rate_constants!(out, gas; do_irreversible=true) == ref
    end

    @testset "kinetics: matrix shapes" begin
        ns = n_species(gas); nr = n_reactions(gas)
        @test size(net_production_rates_ddX(gas)) == (ns, ns)
        @test size(net_rates_of_progress_ddX(gas)) == (nr, ns)
        # Same number of elements, transposed shape: still a size mismatch.
        @test_throws DimensionMismatch net_rates_of_progress_ddX!(zeros(ns, nr), gas)
    end

    @testset "transport: $f" for (f, f!) in (
            (mix_diff_coeffs, mix_diff_coeffs!),
            (thermal_diff_coeffs, thermal_diff_coeffs!),
            (binary_diff_coeffs, binary_diff_coeffs!),
        )
        check_inplace(f, f!, gas)
    end

    @testset "views: Kinetics, ThermoPhase and Transport" begin
        ns = n_species(gas)
        @test net_production_rates!(zeros(ns), kinetics(gas)) == net_production_rates(gas)
        @test mole_fractions!(zeros(ns), thermo(gas)) == mole_fractions(gas)
        @test mix_diff_coeffs!(zeros(ns), transport(gas)) == mix_diff_coeffs(gas)
        @test binary_diff_coeffs!(zeros(ns, ns), transport(gas)) == binary_diff_coeffs(gas)
    end
    close!(gas)

    @testset "transport: multicomponent" begin
        gm = Solution("gri30.yaml"; transport="multicomponent")
        set_TPX!(gm, 1200.0, one_atm, "CH4:1, O2:2, N2:7.52")
        check_inplace(multi_diff_coeffs, multi_diff_coeffs!, gm)
        close!(gm)
    end

    @testset "reactors" begin
        g = Solution("gri30.yaml")
        set_TPX!(g, 1000.0, one_atm, "H2:2, O2:1, N2:4")
        r = IdealGasReactor(g)
        net = ReactorNet(r)
        advance!(net, 1e-4)
        check_inplace(mass_fractions, mass_fractions!, r)
        check_inplace(state, state!, net)
        close!(g)
    end
end

@testset "Setters accept views" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 1500.0, one_atm, "CH4:1, O2:2, N2:7.52")
    X = mole_fractions(gas); Y = mass_fractions(gas)
    M = repeat(Y, 1, 3)
    T, P = temperature(gas), pressure(gas)

    # Contiguous views are passed to the CLib without a copy.
    for (set!, args, get) in (
            (set_mole_fractions!, (view(repeat(X, 1, 2), :, 2),), mole_fractions),
            (set_mass_fractions!, (view(M, :, 2),), mass_fractions),
            (set_TPX!, (T, P, view(repeat(X, 1, 2), :, 2)), mole_fractions),
            (set_TPY!, (T, P, view(M, :, 2)), mass_fractions),
        )
        @test allocations((v, g) -> set!(g, v...), args, gas) == 0
        @test get(gas) ≈ (get === mole_fractions ? X : Y) rtol=1e-12
    end

    # A strided view is copied, and still sets the right state.
    S = zeros(2 * length(Y)); S[1:2:end] .= Y
    set_mass_fractions!(gas, view(S, 1:2:length(S)))
    @test mass_fractions(gas) ≈ Y rtol=1e-12

    # SolutionArray rows are restored from views of its storage.
    sa = SolutionArray(gas, 3)
    @test allocations(restore!, sa, 2) == 0
    close!(gas)
end
