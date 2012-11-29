function v = thermalDiffCoeffs(a)
%THERMALDIFFCOEFFS  Thermal diffusion coefficients.
%
%    dt = thermalDiffCoeffs(gas)
%
%    returns in column vector dt the thermal diffusion
%    coefficients. Object 'gas' must belong to a class derived from
%    Transport, and that was constructed by specifying the 'Multi'
%    option. If 'Multi' was not specified, the returned values will
%    all be zero.
%
v = trans_get(a.id, 12, nSpecies(a.th));
