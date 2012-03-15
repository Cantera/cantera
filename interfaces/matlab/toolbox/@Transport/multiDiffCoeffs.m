function v = multiDiffCoeffs(a)
%MULTIDIFFCOEFFS  Multicomponent diffusion coefficients (m^2/s).
%
%    d = multiDiffCoeffs(gas)
%
%    returns in d the array of multicomponent diffusion
%    coefficients. Object 'gas' must belong to a class derived from
%    Transport, and that was constructed by specifying the 'Multi'
%    option. If 'Multi' was not specified, you will get the error message
%
%    **** Method getMultiDiffCoeffs not implemented. ****
%
%    In this case, try method 'mixDiffCoeffs', or create a new gas
%    mixture model that uses a multicomponent transport manager,
%    for example:
%
%        gas = GRI30('Multi')
%
v = trans_get(a.id, 22, nSpecies(a.th));

