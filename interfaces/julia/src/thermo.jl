# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Thermodynamic state, composition and property accessors.
#
# Every function accepts either a `Solution` or a `ThermoPhase` (collectively
# `ThermoLike`) and operates on the underlying ThermoPhase handle.  Scalar
# getters return `Float64`; array getters return freshly allocated
# `Vector{Float64}` and have in-place `!` variants for hot loops.

# ---- sizes and names --------------------------------------------------------

"Number of chemical species in the phase."
n_species(g::ThermoLike) = Int(check(LibCantera.thermo_nSpecies(_thermo_handle(g))))

"Number of elements in the phase."
n_elements(g::ThermoLike) = Int(check(LibCantera.thermo_nElements(_thermo_handle(g))))

"Name of species `k` (1-based)."
function species_name(g::ThermoLike, k::Integer)
    h = _thermo_handle(g)
    return get_string((n, b) -> LibCantera.thermo_speciesName(h, Int32(k - 1), n, b))
end

"Vector of all species names."
species_names(g::ThermoLike) = [species_name(g, k) for k in 1:n_species(g)]

"""
    species_index(gas, name) -> Int

1-based index of species `name`, or `0` if it is not present.
"""
function species_index(g::ThermoLike, nm::AbstractString)
    idx = LibCantera.thermo_speciesIndex(_thermo_handle(g), nm)
    idx == -1 && return 0            # not found (CLib returns npos-as-(-1))
    return Int(idx) + 1
end

"Name of element `m` (1-based)."
function element_name(g::ThermoLike, m::Integer)
    h = _thermo_handle(g)
    return get_string((n, b) -> LibCantera.thermo_elementName(h, Int32(m - 1), n, b))
end

"Vector of all element names."
element_names(g::ThermoLike) = [element_name(g, m) for m in 1:n_elements(g)]

# ---- scalar state -----------------------------------------------------------

"Temperature [K]."
temperature(g::ThermoLike) = checkd(LibCantera.thermo_temperature(_thermo_handle(g)))

"Pressure [Pa]."
pressure(g::ThermoLike) = checkd(LibCantera.thermo_pressure(_thermo_handle(g)))

"Mass density [kg/m^3]."
density(g::ThermoLike) = checkd(LibCantera.thermo_density(_thermo_handle(g)))

"Molar density [kmol/m^3]."
molar_density(g::ThermoLike) = checkd(LibCantera.thermo_molarDensity(_thermo_handle(g)))

"Mean molecular weight [kg/kmol]."
mean_molecular_weight(g::ThermoLike) =
    checkd(LibCantera.thermo_meanMolecularWeight(_thermo_handle(g)))

# ---- molar / specific properties -------------------------------------------

for (jl, c, doc) in (
        (:enthalpy_mass,        :thermo_enthalpy_mass,  "Specific enthalpy [J/kg]."),
        (:enthalpy_mole,        :thermo_enthalpy_mole,  "Molar enthalpy [J/kmol]."),
        (:internal_energy_mass, :thermo_intEnergy_mass,
            "Specific internal energy [J/kg]."),
        (:internal_energy_mole, :thermo_intEnergy_mole,
            "Molar internal energy [J/kmol]."),
        (:entropy_mass,         :thermo_entropy_mass,   "Specific entropy [J/kg/K]."),
        (:entropy_mole,         :thermo_entropy_mole,   "Molar entropy [J/kmol/K]."),
        (:gibbs_mass,           :thermo_gibbs_mass,     "Specific Gibbs free energy [J/kg]."),
        (:gibbs_mole,           :thermo_gibbs_mole,     "Molar Gibbs free energy [J/kmol]."),
        (:cp_mass,              :thermo_cp_mass,
            "Specific heat capacity at constant pressure [J/kg/K]."),
        (:cp_mole,              :thermo_cp_mole,
            "Molar heat capacity at constant pressure [J/kmol/K]."),
        (:cv_mass,              :thermo_cv_mass,
            "Specific heat capacity at constant volume [J/kg/K]."),
        (:cv_mole,              :thermo_cv_mole,
            "Molar heat capacity at constant volume [J/kmol/K]."),
    )
    @eval begin
        $jl(g::ThermoLike) = checkd(LibCantera.$c(_thermo_handle(g)))
        @doc $doc $jl
    end
end

"Isothermal compressibility [1/Pa]."
isothermal_compressibility(g::ThermoLike) =
    checkd(LibCantera.thermo_isothermalCompressibility(_thermo_handle(g)))

"Thermal (volumetric) expansion coefficient [1/K]."
thermal_expansion_coeff(g::ThermoLike) =
    checkd(LibCantera.thermo_thermalExpansionCoeff(_thermo_handle(g)))

"Maximum temperature [K] for which the phase's thermo data are valid."
max_temp(g::ThermoLike) =
    checkd(LibCantera.thermo_maxTemp(_thermo_handle(g), Int32(-1)))

"Minimum temperature [K] for which the phase's thermo data are valid."
min_temp(g::ThermoLike) =
    checkd(LibCantera.thermo_minTemp(_thermo_handle(g), Int32(-1)))

"Reference pressure [Pa] used for standard-state thermo data."
reference_pressure(g::ThermoLike) =
    checkd(LibCantera.thermo_refPressure(_thermo_handle(g)))

"Electric potential [V] of the phase."
electric_potential(g::ThermoLike) =
    checkd(LibCantera.thermo_electricPotential(_thermo_handle(g)))

"Set the electric potential [V] of the phase."
function set_electric_potential!(g::ThermoLike, v)
    check(LibCantera.thermo_setElectricPotential(_thermo_handle(g), Float64(v)))
    return g
end

"Speed of sound [m/s] (`sqrt(cp/cv · p/ρ)`)."
sound_speed(g::ThermoLike) = sqrt(cp_mass(g) / cv_mass(g) * pressure(g) / density(g))

"Specific volume [m^3/kg]."
volume_mass(g::ThermoLike) = 1.0 / density(g)

"Molar volume [m^3/kmol]."
volume_mole(g::ThermoLike) = mean_molecular_weight(g) / density(g)

# ---- critical / saturation properties (raise for ideal-gas phases) ----------

"Critical temperature [K]."
critical_temperature(g::ThermoLike) =
    checkd(LibCantera.thermo_critTemperature(_thermo_handle(g)))
"Critical pressure [Pa]."
critical_pressure(g::ThermoLike) =
    checkd(LibCantera.thermo_critPressure(_thermo_handle(g)))
"Critical density [kg/m^3]."
critical_density(g::ThermoLike) =
    checkd(LibCantera.thermo_critDensity(_thermo_handle(g)))
"Vapor fraction (quality) of a two-phase state."
vapor_fraction(g::ThermoLike) =
    checkd(LibCantera.thermo_vaporFraction(_thermo_handle(g)))
"Saturation temperature [K] at pressure `p` [Pa]."
sat_temperature(g::ThermoLike, p) =
    checkd(LibCantera.thermo_satTemperature(_thermo_handle(g), Float64(p)))
"Saturation pressure [Pa] at temperature `T` [K]."
sat_pressure(g::ThermoLike, T) =
    checkd(LibCantera.thermo_satPressure(_thermo_handle(g), Float64(T)))

# ---- elements / atoms -------------------------------------------------------

"""
    element_index(gas, name) -> Int

1-based index of element `name`, or `0` if it is not present.
"""
function element_index(g::ThermoLike, nm::AbstractString)
    idx = LibCantera.thermo_elementIndex(_thermo_handle(g), nm)
    idx == -1 && return 0
    return Int(idx) + 1
end

"Number of atoms of element `m` (1-based) in species `k` (1-based)."
n_atoms(g::ThermoLike, k::Integer, m::Integer) =
    checkd(LibCantera.thermo_nAtoms(_thermo_handle(g), Int32(k - 1), Int32(m - 1)))

"Atomic weights of all elements [kg/kmol]."
atomic_weights(g::ThermoLike) =
    get_array(n_elements(g),
        (n, b) -> LibCantera.thermo_atomicWeights(_thermo_handle(g), n, b))

"Atomic weight of element `m` (1-based) [kg/kmol]."
atomic_weight(g::ThermoLike, m::Integer) = atomic_weights(g)[m]

"Species electric charges (per elementary charge)."
charges(g::ThermoLike) =
    get_array(n_species(g),
        (n, b) -> LibCantera.thermo_getCharges(_thermo_handle(g), n, b))

"""
    elemental_mole_fraction(gas, m) -> Float64

Mole fraction of element `m` (1-based) — moles of that element per mole of all
elements in the mixture.
"""
function elemental_mole_fraction(g::ThermoLike, m::Integer)
    X = mole_fractions(g); ne = n_elements(g); nsp = n_species(g)
    num = sum(n_atoms(g, k, m) * X[k] for k in 1:nsp)
    den = sum(sum(n_atoms(g, k, j) * X[k] for k in 1:nsp) for j in 1:ne)
    return num / den
end

"""
    elemental_mass_fraction(gas, m) -> Float64

Mass fraction of element `m` (1-based) in the mixture.
"""
function elemental_mass_fraction(g::ThermoLike, m::Integer)
    X = mole_fractions(g); nsp = n_species(g)
    num = atomic_weight(g, m) * sum(n_atoms(g, k, m) * X[k] for k in 1:nsp)
    return num / mean_molecular_weight(g)
end

"""
    equivalence_ratio(gas) -> Float64

Fuel/air equivalence ratio from the elemental composition (C, H, O, S), using
complete-combustion available-oxygen accounting: `φ = (2·Z_C + Z_H/2 + 2·Z_S) / Z_O`.
"""
function equivalence_ratio(g::ThermoLike)
    Z(sym) = (i = element_index(g, sym); i == 0 ? 0.0 : elemental_mole_fraction(g, i))
    ZO = Z("O")
    ZO == 0 && return Inf
    return (2Z("C") + Z("H") / 2 + 2Z("S")) / ZO
end

# ---- composition setters ----------------------------------------------------

"""
    set_equivalence_ratio!(gas, phi, fuel, oxidizer)

Set the mixture composition to the given equivalence ratio `phi` for the given
`fuel` and `oxidizer` compositions (name-value strings), holding T and p.
"""
function set_equivalence_ratio!(g::ThermoLike, phi, fuel::AbstractString,
                                oxidizer::AbstractString)
    check(LibCantera.thermo_setEquivalenceRatio(_thermo_handle(g), Float64(phi),
                                                fuel, oxidizer))
    return g
end

# ---- standard-state (per species) -------------------------------------------

for (jl, c, doc) in (
        (:standard_enthalpies_RT,   :thermo_getEnthalpy_RT,  "enthalpies h°_k/RT"),
        (:standard_entropies_R,     :thermo_getEntropy_R,    "entropies s°_k/R"),
        (:standard_gibbs_RT,        :thermo_getGibbs_RT,     "Gibbs energies g°_k/RT"),
        (:standard_int_energies_RT, :thermo_getIntEnergy_RT,
         "internal energies u°_k/RT"),
        (:standard_cp_R,            :thermo_getCp_R,         "heat capacities cp°_k/R"),
    )
    @eval begin
        @doc $("Nondimensional standard-state species " * doc * ".")
        $jl(g::ThermoLike) =
            get_array(n_species(g), (n, b) -> LibCantera.$c(_thermo_handle(g), n, b))
    end
end

# ---- composition ------------------------------------------------------------

"Molecular weights of all species [kg/kmol]."
molecular_weights(g::ThermoLike) =
    get_array(n_species(g),
        (n, b) -> LibCantera.thermo_getMolecularWeights(_thermo_handle(g), n, b))

"In-place [`molecular_weights`](@ref)."
molecular_weights!(out, g::ThermoLike) =
    get_array!(out,
        (n, b) -> LibCantera.thermo_getMolecularWeights(_thermo_handle(g), n, b))

"Mole fractions of all species."
mole_fractions(g::ThermoLike) =
    get_array(n_species(g),
        (n, b) -> LibCantera.thermo_getMoleFractions(_thermo_handle(g), n, b))

"In-place [`mole_fractions`](@ref)."
mole_fractions!(out, g::ThermoLike) =
    get_array!(out,
        (n, b) -> LibCantera.thermo_getMoleFractions(_thermo_handle(g), n, b))

"Mass fractions of all species."
mass_fractions(g::ThermoLike) =
    get_array(n_species(g),
        (n, b) -> LibCantera.thermo_getMassFractions(_thermo_handle(g), n, b))

"In-place [`mass_fractions`](@ref)."
mass_fractions!(out, g::ThermoLike) =
    get_array!(out,
        (n, b) -> LibCantera.thermo_getMassFractions(_thermo_handle(g), n, b))

"Species concentrations [kmol/m^3]."
concentrations(g::ThermoLike) =
    get_array(n_species(g),
        (n, b) -> LibCantera.thermo_getConcentrations(_thermo_handle(g), n, b))

# ---- partial molar / potentials --------------------------------------------

for (jl, c, doc) in (
        (:partial_molar_enthalpies, :thermo_getPartialMolarEnthalpies,
            "Partial molar enthalpies of the species [J/kmol]."),
        (:partial_molar_entropies,  :thermo_getPartialMolarEntropies,
            "Partial molar entropies of the species [J/kmol/K]."),
        (:partial_molar_int_energies, :thermo_getPartialMolarIntEnergies,
            "Partial molar internal energies of the species [J/kmol]."),
        (:partial_molar_cp,         :thermo_getPartialMolarCp,
            "Partial molar heat capacities at constant pressure [J/kmol/K]."),
        (:partial_molar_volumes,    :thermo_getPartialMolarVolumes,
            "Partial molar volumes of the species [m^3/kmol]."),
        (:chemical_potentials,      :thermo_getChemPotentials,
            "Chemical potentials of the species [J/kmol]."),
        (:electrochemical_potentials, :thermo_getElectrochemPotentials,
            "Electrochemical potentials of the species [J/kmol]."),
    )
    bang = Symbol(jl, :!)
    bang_doc = "In-place [`$jl`](@ref)."
    @eval begin
        $jl(g::ThermoLike) =
            get_array(n_species(g), (n, b) -> LibCantera.$c(_thermo_handle(g), n, b))
        @doc $doc $jl
        $bang(out, g::ThermoLike) =
            get_array!(out, (n, b) -> LibCantera.$c(_thermo_handle(g), n, b))
        @doc $bang_doc $bang
    end
end

# ---- composition setters ----------------------------------------------------

"""
    set_mole_fractions!(gas, X)

Set mole fractions from a numeric vector or a composition string
(e.g. `"CH4:1, O2:2"`).  The values are normalized by Cantera.
"""
function set_mole_fractions!(g::ThermoLike, X::AbstractVector{<:Real})
    x = as_f64(X)
    check(LibCantera.thermo_setMoleFractions(_thermo_handle(g), Int32(length(x)),
                                             pointer(x)))
    return g
end
function set_mole_fractions!(g::ThermoLike, X::AbstractString)
    check(LibCantera.thermo_setMoleFractionsByName(_thermo_handle(g), X))
    return g
end

"""
    set_mass_fractions!(gas, Y)

Set mass fractions from a numeric vector or a composition string.
"""
function set_mass_fractions!(g::ThermoLike, Y::AbstractVector{<:Real})
    y = as_f64(Y)
    check(LibCantera.thermo_setMassFractions(_thermo_handle(g), Int32(length(y)),
                                             pointer(y)))
    return g
end
function set_mass_fractions!(g::ThermoLike, Y::AbstractString)
    check(LibCantera.thermo_setMassFractionsByName(_thermo_handle(g), Y))
    return g
end

# ---- state setters ----------------------------------------------------------

"Set temperature [K] and pressure [Pa]."
function set_TP!(g::ThermoLike, T, P)
    check(LibCantera.thermo_setState_TP(_thermo_handle(g), Float64(T), Float64(P)))
    return g
end

"""
    set_TPX!(gas, T, p, X)

Set temperature, pressure and mole fractions in one call.  `X` may be a numeric
vector or a composition string.
"""
function set_TPX!(g::ThermoLike, T, P, X::AbstractVector{<:Real})
    x = as_f64(X)
    check(LibCantera.thermo_setState_TPX(_thermo_handle(g), Float64(T), Float64(P),
                                         Int32(length(x)), pointer(x)))
    return g
end
function set_TPX!(g::ThermoLike, T, P, X::AbstractString)
    check(LibCantera.thermo_setState_TPX_byName(_thermo_handle(g), Float64(T),
                                                 Float64(P), X))
    return g
end

"""
    set_TPY!(gas, T, p, Y)

Set temperature, pressure and mass fractions in one call.  `Y` may be a numeric
vector or a composition string.
"""
function set_TPY!(g::ThermoLike, T, P, Y::AbstractVector{<:Real})
    y = as_f64(Y)
    check(LibCantera.thermo_setState_TPY(_thermo_handle(g), Float64(T), Float64(P),
                                         Int32(length(y)), pointer(y)))
    return g
end
function set_TPY!(g::ThermoLike, T, P, Y::AbstractString)
    check(LibCantera.thermo_setState_TPY_byName(_thermo_handle(g), Float64(T),
                                                 Float64(P), Y))
    return g
end

# Two-property setters that map directly onto CLib `setState_XY`.
for (jl, c, a, b, doc) in (
        (:set_HP!, :thermo_setState_HP, :h, :p,
            "Set specific enthalpy `h` [J/kg] and pressure `p` [Pa]."),
        (:set_UV!, :thermo_setState_UV, :u, :v,
            "Set specific internal energy `u` [J/kg] and specific volume `v` [m^3/kg]."),
        (:set_SP!, :thermo_setState_SP, :s, :p,
            "Set specific entropy `s` [J/kg/K] and pressure `p` [Pa]."),
        (:set_SV!, :thermo_setState_SV, :s, :v,
            "Set specific entropy `s` [J/kg/K] and specific volume `v` [m^3/kg]."),
        (:set_TD!, :thermo_setState_TD, :T, :rho,
            "Set temperature `T` [K] and density `rho` [kg/m^3]."),
        (:set_DP!, :thermo_setState_DP, :rho, :p,
            "Set density `rho` [kg/m^3] and pressure `p` [Pa]."),
    )
    @eval begin
        function $jl(g::ThermoLike, $a, $b)
            check(LibCantera.$c(_thermo_handle(g), Float64($a), Float64($b)))
            return g
        end
        @doc $doc $jl
    end
end

"""
    equilibrate!(gas, XY; solver="auto", rtol=1e-9, max_steps=1000,
                 max_iter=100, estimate_equil=0)

Set the phase to a state of chemical equilibrium holding the property pair `XY`
(e.g. `"TP"`, `"HP"`, `"UV"`) fixed.
"""
function equilibrate!(g::ThermoLike, XY::AbstractString; solver="auto",
                      rtol=1e-9, max_steps=1000, max_iter=100, estimate_equil=0)
    check(LibCantera.thermo_equilibrate(_thermo_handle(g), XY, solver, Float64(rtol),
                                        Int32(max_steps), Int32(max_iter),
                                        Int32(estimate_equil)))
    return g
end

"""
    report(gas; show_thermo=true, threshold=1e-14) -> String

Return Cantera's formatted state report for the phase.
"""
function report(g::ThermoLike; show_thermo::Bool=true, threshold=1e-14)
    h = _thermo_handle(g)
    return get_string((n, b) -> LibCantera.thermo_report(h, Int32(show_thermo),
                                                         Float64(threshold), n, b))
end

export report

export n_species, n_elements, species_name, species_names, species_index,
       element_name, element_names,
       mole_fractions, mole_fractions!, mass_fractions, mass_fractions!,
       set_mole_fractions!, set_mass_fractions!, molecular_weights,
       molecular_weights!, concentrations, mean_molecular_weight

export temperature, pressure, density, molar_density,
       enthalpy_mass, enthalpy_mole, internal_energy_mass, internal_energy_mole,
       entropy_mass, entropy_mole, gibbs_mass, gibbs_mole,
       cp_mass, cp_mole, cv_mass, cv_mole,
       isothermal_compressibility, thermal_expansion_coeff,
       reference_pressure, electric_potential, set_electric_potential!,
       standard_enthalpies_RT, standard_entropies_R, standard_gibbs_RT,
       standard_int_energies_RT, standard_cp_R,
       sound_speed, volume_mass, volume_mole,
       critical_temperature, critical_pressure, critical_density,
       vapor_fraction, sat_temperature, sat_pressure,
       element_index, atomic_weights, atomic_weight, charges,
       elemental_mole_fraction, elemental_mass_fraction, equivalence_ratio,
       partial_molar_enthalpies, partial_molar_entropies,
       partial_molar_int_energies, partial_molar_cp, partial_molar_volumes,
       chemical_potentials, electrochemical_potentials,
       partial_molar_enthalpies!, partial_molar_entropies!, partial_molar_cp!,
       chemical_potentials!

export set_TP!, set_TPX!, set_TPY!, set_HP!, set_UV!, set_SP!, set_SV!,
       set_TD!, set_DP!, equilibrate!, set_equivalence_ratio!
