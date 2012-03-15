function v = mixDiffCoeffs(a)
%MIXDIFFCOEFFS  Mixture-averaged diffusion coefficients (m^2/s).
%
%    d = mixDiffCoeffs(gas)
%
%    returns in column vector d the mixture-averaged diffusion
%    coefficients. Object 'gas' must belong to a class derived from
%    Transport, and that was constructed by specifying the 'Mix'
%    option. If 'Mix' was not specified, you will get the error message
%
%    **** Method getMixDiffCoeffs not implemented. ****
%
%    In this case, try method 'multiDiffCoeffs', or create a new gas
%    mixture model that uses a mixture-averaged transport manager,
%    for example:
%
%        gas = GRI30('Mix')
%
v = trans_get(a.id, 11, nSpecies(a.th));
