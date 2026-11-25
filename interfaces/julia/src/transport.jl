# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Transport-property accessors.  Accept a `Solution` or a `Transport`.

"Dynamic (mixture) viscosity [Pa*s]."
viscosity(g::TransportLike) = checkd(LibCantera.trans_viscosity(_transport_handle(g)))

"Mixture thermal conductivity [W/m/K]."
thermal_conductivity(g::TransportLike) =
    checkd(LibCantera.trans_thermalConductivity(_transport_handle(g)))

"Electrical conductivity [S/m]."
electrical_conductivity(g::TransportLike) =
    checkd(LibCantera.trans_electricalConductivity(_transport_handle(g)))

# Species count for the per-species transport arrays.  A `Transport` view has no
# thermo handle of its own, so it is taken from the owning `Solution`.
_transport_nsp(g::Solution) = n_species(g)
_transport_nsp(t::Transport) = n_species(t.parent)

"Mixture-averaged diffusion coefficients [m^2/s], one per species."
mix_diff_coeffs(g::Solution) =
    get_array(n_species(g),
              (n, b) -> LibCantera.trans_getMixDiffCoeffs(_transport_handle(g), n, b))

mix_diff_coeffs!(out, g::TransportLike) =
    get_array!(out, _transport_nsp(g),
               (n, b) -> LibCantera.trans_getMixDiffCoeffs(_transport_handle(g), n, b))
@doc inplace_doc(:mix_diff_coeffs) mix_diff_coeffs!

"Thermal diffusion coefficients [kg/m/s], one per species."
thermal_diff_coeffs(g::Solution) =
    get_array(n_species(g),
              (n, b) -> LibCantera.trans_getThermalDiffCoeffs(
                  _transport_handle(g), n, b))

thermal_diff_coeffs!(out, g::TransportLike) =
    get_array!(out, _transport_nsp(g),
               (n, b) -> LibCantera.trans_getThermalDiffCoeffs(
                   _transport_handle(g), n, b))
@doc inplace_doc(:thermal_diff_coeffs) thermal_diff_coeffs!

"""
    binary_diff_coeffs(gas) -> Matrix{Float64}

`nsp x nsp` matrix of binary diffusion coefficients [m^2/s].
"""
binary_diff_coeffs(g::Solution) =
    binary_diff_coeffs!(Matrix{Float64}(undef, n_species(g), n_species(g)), g)

function binary_diff_coeffs!(out, g::TransportLike)
    nsp = _transport_nsp(g)
    return get_matrix!(out, nsp, nsp,
        (len, b) -> LibCantera.trans_getBinaryDiffCoeffs(
            _transport_handle(g), Int32(nsp), len, b))
end
@doc inplace_doc(:binary_diff_coeffs) binary_diff_coeffs!

"""
    multi_diff_coeffs(gas) -> Matrix{Float64}

`nsp x nsp` matrix of multicomponent diffusion coefficients [m^2/s].
"""
multi_diff_coeffs(g::Solution) =
    multi_diff_coeffs!(Matrix{Float64}(undef, n_species(g), n_species(g)), g)

function multi_diff_coeffs!(out, g::TransportLike)
    nsp = _transport_nsp(g)
    return get_matrix!(out, nsp, nsp,
        (len, b) -> LibCantera.trans_getMultiDiffCoeffs(
            _transport_handle(g), Int32(nsp), len, b))
end
@doc inplace_doc(:multi_diff_coeffs) multi_diff_coeffs!

export viscosity, thermal_conductivity, electrical_conductivity,
       mix_diff_coeffs, mix_diff_coeffs!, thermal_diff_coeffs, thermal_diff_coeffs!,
       binary_diff_coeffs, binary_diff_coeffs!, multi_diff_coeffs, multi_diff_coeffs!
