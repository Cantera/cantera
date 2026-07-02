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

for (jl, c) in (
        (:enthalpy_mass,        :thermo_enthalpy_mass),
        (:enthalpy_mole,        :thermo_enthalpy_mole),
        (:internal_energy_mass, :thermo_intEnergy_mass),
        (:internal_energy_mole, :thermo_intEnergy_mole),
        (:entropy_mass,         :thermo_entropy_mass),
        (:entropy_mole,         :thermo_entropy_mole),
        (:gibbs_mass,           :thermo_gibbs_mass),
        (:gibbs_mole,           :thermo_gibbs_mole),
        (:cp_mass,              :thermo_cp_mass),
        (:cp_mole,              :thermo_cp_mole),
        (:cv_mass,              :thermo_cv_mass),
        (:cv_mole,              :thermo_cv_mole),
    )
    @eval $jl(g::ThermoLike) = checkd(LibCantera.$c(_thermo_handle(g)))
end

"Isothermal compressibility [1/Pa]."
isothermal_compressibility(g::ThermoLike) =
    checkd(LibCantera.thermo_isothermalCompressibility(_thermo_handle(g)))

"Thermal (volumetric) expansion coefficient [1/K]."
thermal_expansion_coeff(g::ThermoLike) =
    checkd(LibCantera.thermo_thermalExpansionCoeff(_thermo_handle(g)))

# The CLib getters take a species index; passing -1 (C++ `npos`) requests the
# phase-wide limit, matching the default `maxTemp()`/`minTemp()` in Python.
"Maximum temperature [K] for which the phase's thermo data are valid."
max_temp(g::ThermoLike) = checkd(LibCantera.thermo_maxTemp(_thermo_handle(g), Int32(-1)))

"Minimum temperature [K] for which the phase's thermo data are valid."
min_temp(g::ThermoLike) = checkd(LibCantera.thermo_minTemp(_thermo_handle(g), Int32(-1)))

# ---- composition ------------------------------------------------------------

"Molecular weights of all species [kg/kmol]."
molecular_weights(g::ThermoLike) =
    get_array(n_species(g), (n, b) -> LibCantera.thermo_getMolecularWeights(_thermo_handle(g), n, b))

"In-place [`molecular_weights`](@ref)."
molecular_weights!(out, g::ThermoLike) =
    get_array!(out, (n, b) -> LibCantera.thermo_getMolecularWeights(_thermo_handle(g), n, b))

"Mole fractions of all species."
mole_fractions(g::ThermoLike) =
    get_array(n_species(g), (n, b) -> LibCantera.thermo_getMoleFractions(_thermo_handle(g), n, b))

"In-place [`mole_fractions`](@ref)."
mole_fractions!(out, g::ThermoLike) =
    get_array!(out, (n, b) -> LibCantera.thermo_getMoleFractions(_thermo_handle(g), n, b))

"Mass fractions of all species."
mass_fractions(g::ThermoLike) =
    get_array(n_species(g), (n, b) -> LibCantera.thermo_getMassFractions(_thermo_handle(g), n, b))

"In-place [`mass_fractions`](@ref)."
mass_fractions!(out, g::ThermoLike) =
    get_array!(out, (n, b) -> LibCantera.thermo_getMassFractions(_thermo_handle(g), n, b))

"Species concentrations [kmol/m^3]."
concentrations(g::ThermoLike) =
    get_array(n_species(g), (n, b) -> LibCantera.thermo_getConcentrations(_thermo_handle(g), n, b))

# ---- partial molar / potentials --------------------------------------------

for (jl, c) in (
        (:partial_molar_enthalpies, :thermo_getPartialMolarEnthalpies),
        (:partial_molar_entropies,  :thermo_getPartialMolarEntropies),
        (:partial_molar_int_energies, :thermo_getPartialMolarIntEnergies),
        (:partial_molar_cp,         :thermo_getPartialMolarCp),
        (:partial_molar_volumes,    :thermo_getPartialMolarVolumes),
        (:chemical_potentials,      :thermo_getChemPotentials),
        (:electrochemical_potentials, :thermo_getElectrochemPotentials),
    )
    bang = Symbol(jl, :!)
    @eval begin
        $jl(g::ThermoLike) =
            get_array(n_species(g), (n, b) -> LibCantera.$c(_thermo_handle(g), n, b))
        $bang(out, g::ThermoLike) =
            get_array!(out, (n, b) -> LibCantera.$c(_thermo_handle(g), n, b))
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
    check(LibCantera.thermo_setMoleFractions(_thermo_handle(g), Int32(length(x)), pointer(x)))
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
    check(LibCantera.thermo_setMassFractions(_thermo_handle(g), Int32(length(y)), pointer(y)))
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
    check(LibCantera.thermo_setState_TPX_byName(_thermo_handle(g), Float64(T), Float64(P), X))
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
    check(LibCantera.thermo_setState_TPY_byName(_thermo_handle(g), Float64(T), Float64(P), Y))
    return g
end

# Two-property setters that map directly onto CLib `setState_XY`.
for (jl, c, a, b) in (
        (:set_HP!, :thermo_setState_HP, :h, :p),
        (:set_UV!, :thermo_setState_UV, :u, :v),
        (:set_SP!, :thermo_setState_SP, :s, :p),
        (:set_SV!, :thermo_setState_SV, :s, :v),
        (:set_TD!, :thermo_setState_TD, :T, :rho),
        (:set_DP!, :thermo_setState_DP, :rho, :p),
    )
    @eval function $jl(g::ThermoLike, $a, $b)
        check(LibCantera.$c(_thermo_handle(g), Float64($a), Float64($b)))
        return g
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
