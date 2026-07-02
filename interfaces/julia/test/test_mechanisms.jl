using Cantera
using Test

# A second mechanism (h2o2.yaml, 10 species / 29 reactions) to complement the
# gri30-based reference checks and guard against mechanism-specific assumptions.
# Reference values from Python Cantera 3.x at 1000 K, 1 atm, H2:2/O2:1/N2:4.

@testset "h2o2 mechanism vs reference" begin
    gas = Solution("h2o2.yaml")
    set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")

    @test n_species(gas) == 10
    @test n_reactions(gas) == 29
    @test density(gas) ≈ 0.25780918724920693 rtol=1e-8
    @test cp_mass(gas) ≈ 1527.8760405617072 rtol=1e-8
    @test enthalpy_mass(gas) ≈ 1012650.3485874726 rtol=1e-8
    @test viscosity(gas) ≈ 4.199970957871177e-5 rtol=1e-6
    @test thermal_conductivity(gas) ≈ 0.1317893909081434 rtol=1e-6

    @test species_index(gas, "OH") == 5           # 1-based (Python 0-based 4)
    @test reaction_equations(gas)[1] == "2 O + M <=> O2 + M"

    # consistency identities that must hold for any mechanism
    @test net_production_rates(gas) ≈ creation_rates(gas) .- destruction_rates(gas) rtol=1e-8
    @test net_rates_of_progress(gas) ≈
          forward_rates_of_progress(gas) .- reverse_rates_of_progress(gas) rtol=1e-10
    close!(gas)
end
