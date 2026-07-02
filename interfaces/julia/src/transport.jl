# Transport-property accessors.  Accept a `Solution` or a `Transport`.

"Dynamic (mixture) viscosity [Pa*s]."
viscosity(g::TransportLike) = checkd(LibCantera.trans_viscosity(_transport_handle(g)))

"Mixture thermal conductivity [W/m/K]."
thermal_conductivity(g::TransportLike) =
    checkd(LibCantera.trans_thermalConductivity(_transport_handle(g)))

"Electrical conductivity [S/m]."
electrical_conductivity(g::TransportLike) =
    checkd(LibCantera.trans_electricalConductivity(_transport_handle(g)))

"Mixture-averaged diffusion coefficients [m^2/s], one per species."
mix_diff_coeffs(g::Solution) =
    get_array(n_species(g),
              (n, b) -> LibCantera.trans_getMixDiffCoeffs(_transport_handle(g), n, b))

"In-place [`mix_diff_coeffs`](@ref)."
mix_diff_coeffs!(out, g::TransportLike) =
    get_array!(out, (n, b) -> LibCantera.trans_getMixDiffCoeffs(_transport_handle(g), n, b))

"Thermal diffusion coefficients [kg/m/s], one per species."
thermal_diff_coeffs(g::Solution) =
    get_array(n_species(g),
              (n, b) -> LibCantera.trans_getThermalDiffCoeffs(_transport_handle(g), n, b))

"""
    binary_diff_coeffs(gas) -> Matrix{Float64}

`nsp x nsp` matrix of binary diffusion coefficients [m^2/s].
"""
function binary_diff_coeffs(g::Solution)
    nsp = n_species(g)
    buf = Vector{Float64}(undef, nsp * nsp)
    check(LibCantera.trans_getBinaryDiffCoeffs(_transport_handle(g), Int32(nsp),
                                               Int32(length(buf)), pointer(buf)))
    return reshape(buf, nsp, nsp)
end

"""
    multi_diff_coeffs(gas) -> Matrix{Float64}

`nsp x nsp` matrix of multicomponent diffusion coefficients [m^2/s].
"""
function multi_diff_coeffs(g::Solution)
    nsp = n_species(g)
    buf = Vector{Float64}(undef, nsp * nsp)
    check(LibCantera.trans_getMultiDiffCoeffs(_transport_handle(g), Int32(nsp),
                                              Int32(length(buf)), pointer(buf)))
    return reshape(buf, nsp, nsp)
end
