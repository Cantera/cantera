"""
    Cantera

Julia interface to [Cantera](https://cantera.org), built on top of
Cantera's generated CLib API.

```julia
using Cantera
gas = Solution("gri30.yaml")
set_TPY!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
temperature(gas), pressure(gas), net_production_rates(gas)
```
"""
module Cantera

include("LibCantera.jl")
using .LibCantera: LibCantera, libcantera

include("errors.jl")
include("handles.jl")
include("solution.jl")
include("thermo.jl")
include("kinetics.jl")
include("transport.jl")
include("reaction.jl")
include("func1.jl")
include("reactor.jl")
include("reactornet.jl")
include("connectors.jl")
include("onedim.jl")
include("multiphase.jl")
include("rdiag.jl")
include("solutionarray.jl")
include("utils.jl")

# Register Cantera's bundled data directory so mechanisms can be loaded by name.
function __init__()
    d = LibCantera.default_data_directory()
    if d !== nothing
        try
            LibCantera.ct_addDataDirectory(d)
        catch
        end
    end
    return nothing
end

# ---- public API -------------------------------------------------------------
export CanteraObject, CanteraError, close!
export Solution, ThermoPhase, Kinetics, Transport, Reaction
export thermo, kinetics, transport, name, transport_model, report

# constants
export one_atm, gas_constant, avogadro
export cantera_version, git_commit, add_data_directory, suppress_thermo_warnings,
       reset_storage

# thermo: sizes / names / composition
export n_species, n_elements, species_name, species_names, species_index,
       element_name, element_names,
       mole_fractions, mole_fractions!, mass_fractions, mass_fractions!,
       set_mole_fractions!, set_mass_fractions!, molecular_weights,
       molecular_weights!, concentrations, mean_molecular_weight

# thermo: scalar state and properties
export temperature, pressure, density, molar_density,
       enthalpy_mass, enthalpy_mole, internal_energy_mass, internal_energy_mole,
       entropy_mass, entropy_mole, gibbs_mass, gibbs_mole,
       cp_mass, cp_mole, cv_mass, cv_mole,
       isothermal_compressibility, thermal_expansion_coeff,
       partial_molar_enthalpies, partial_molar_entropies,
       partial_molar_int_energies, partial_molar_cp, partial_molar_volumes,
       chemical_potentials, electrochemical_potentials,
       partial_molar_enthalpies!, partial_molar_entropies!, partial_molar_cp!,
       chemical_potentials!

# thermo: state setters
export set_TP!, set_TPX!, set_TPY!, set_HP!, set_UV!, set_SP!, set_SV!,
       set_TD!, set_DP!, equilibrate!

# kinetics
export n_reactions, n_total_species, reaction_equation, reaction_equations,
       is_reversible, reaction, equation, reaction_type, uses_third_body,
       forward_rate_constants, reverse_rate_constants, equilibrium_constants,
       forward_rates_of_progress, reverse_rates_of_progress,
       net_rates_of_progress, net_production_rates, net_production_rates!,
       creation_rates, creation_rates!, destruction_rates, destruction_rates!,
       forward_rates_of_progress!, reverse_rates_of_progress!,
       net_rates_of_progress!, forward_rate_constants!, equilibrium_constants!,
       delta_enthalpy, delta_gibbs, delta_entropy,
       net_production_rates_ddT, net_production_rates_ddP, net_production_rates_ddC,
       creation_rates_ddT, creation_rates_ddP, creation_rates_ddC,
       destruction_rates_ddT, destruction_rates_ddP, destruction_rates_ddC,
       net_rates_of_progress_ddT, net_rates_of_progress_ddP, net_rates_of_progress_ddC,
       net_production_rates_ddX, net_rates_of_progress_ddX,
       net_production_rates_ddT!, net_production_rates_ddP!, net_production_rates_ddC!,
       creation_rates_ddT!, creation_rates_ddP!, creation_rates_ddC!,
       destruction_rates_ddT!, destruction_rates_ddP!, destruction_rates_ddC!,
       net_rates_of_progress_ddT!, net_rates_of_progress_ddP!, net_rates_of_progress_ddC!,
       heat_release_rate

# transport
export viscosity, thermal_conductivity, electrical_conductivity,
       mix_diff_coeffs, mix_diff_coeffs!, thermal_diff_coeffs,
       binary_diff_coeffs, multi_diff_coeffs

# reactors
export Reactor, IdealGasReactor, ConstPressureReactor,
       IdealGasConstPressureReactor, ReactorNet,
       advance!, step!, set_initial_time!, set_max_time_step!,
       set_tolerances!, set_energy_enabled!, set_chemistry_enabled!,
       set_initial_volume!, mass, volume

# reactor connectors, reservoirs, surfaces, sensitivities
export Connector, FlowDevice, Wall, MassFlowController, Valve, PressureController,
       Reservoir, ReactorSurface,
       mass_flow_rate, set_mass_flow_rate!, device_coefficient,
       set_device_coefficient!, set_pressure_function!, set_time_function!,
       set_primary!, expansion_rate, heat_rate, area, set_area!,
       set_heat_transfer_coeff!, set_thermal_resistance!,
       set_expansion_rate_coeff!, set_emissivity!, set_heat_flux!, set_velocity!,
       connector_type, set_name!, add_sensitivity_reaction!, n_sens_params,
       rtol, atol, set_sensitivity_tolerances!, sensitivity, state

# Func1 function objects
export Func1, evaluate, derivative, func_type, constant_function

# MultiPhase mixtures (note: `cp` is available as `Cantera.cp` to avoid the
# clash with `Base.cp`)
export MultiPhase, add_phase!, init!, update_phases!, n_phases,
       set_temperature!, set_pressure!, min_temp, max_temp, charge,
       enthalpy, entropy, gibbs, volume, phase_moles, set_phase_moles!,
       species_moles, element_moles, mole_fraction, n_atoms, set_moles!

# 1-D flames
export Domain1D, FreeFlame, BurnerFlame, set_burner!, burner_mdot,
       n_points, n_components, component_name,
       component_names, domain_type, grid, value, solution_profile, set_profile!,
       set_flat_profile!, setup_uniform_grid!, setup_grid!, set_refine_criteria!,
       solve!, flame_speed, flame_T, flame_X, flame_velocity,
       set_fixed_temperature!, set_inlet!

# reaction-path diagrams
export ReactionPathDiagram, build!, get_dot, get_data, get_log,
       threshold, set_threshold!, bold_threshold, set_bold_threshold!,
       normal_threshold, set_normal_threshold!, label_threshold,
       set_label_threshold!, scale, set_scale!, arrow_width, set_arrow_width!,
       flow_type, set_flow_type!, title, set_title!, font, set_font!,
       show_details, set_show_details!, bold_color, set_bold_color!,
       normal_color, set_normal_color!, dashed_color, set_dashed_color!,
       dot_options, set_dot_options!, display_only!

# SolutionArray (batch states)
export SolutionArray, restore!, set_state!, extract, write_csv

end # module Cantera
