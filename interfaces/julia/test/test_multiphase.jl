using Cantera
using Test

# Reference values cross-checked against Python `cantera.Mixture` (v3.2.0) for a
# single gri30 gas phase at 1200 K, 1 atm, CH4:1 O2:2 N2:7.52.

const CT = Cantera
const MECH = "gri30.yaml"

@testset "MultiPhase" begin
    gas = Solution(MECH)
    set_TPX!(gas, 1200.0, one_atm, "CH4:1, O2:2, N2:7.52")

    mix = CT.MultiPhase([gas => 1.0])
    @test CT.n_phases(mix) == 1
    @test CT.n_species(mix) == 53
    @test CT.n_elements(mix) == 5

    CT.set_temperature!(mix, 1200.0)
    CT.set_pressure!(mix, one_atm)
    @test CT.temperature(mix) ≈ 1200.0
    @test CT.pressure(mix) ≈ one_atm rtol = 1e-10
    @test CT.phase_moles(mix, 1) ≈ 1.0 rtol = 1e-8

    ci = CT.element_index(mix, "C")
    @test ci > 0
    @test CT.element_moles(mix, ci) ≈ 0.09505703422053233 rtol = 1e-6

    CT.equilibrate!(mix, "TP")
    @test CT.temperature(mix) ≈ 1200.0 rtol = 1e-8
    @test CT.phase_moles(mix, 1) ≈ 1.0000017646000376 rtol = 1e-6

    # gibbs == sum(mu_k * moles_k)
    mu = CT.chemical_potentials(mix)
    sm = [CT.species_moles(mix, k) for k in 1:CT.n_species(mix)]
    @test CT.gibbs(mix) ≈ sum(mu .* sm) rtol = 1e-8
    @test CT.gibbs(mix) ≈ -347848546.95592207 rtol = 1e-5

    @test isfinite(CT.enthalpy(mix))
    @test isfinite(CT.entropy(mix))
    @test isfinite(CT.cp(mix))
    @test isfinite(CT.volume(mix))
    @test CT.min_temp(mix) < CT.max_temp(mix)

    @test occursin("phases", sprint(show, mix))

    CT.close!(mix)
    close!(gas)
end
