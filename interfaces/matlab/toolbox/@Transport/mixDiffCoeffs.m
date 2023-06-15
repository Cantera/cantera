function v = mixDiffCoeffs(a)
% MIXDIFFCOEFFS  Get the mixture-averaged diffusion coefficients.
% v = mixDiffCoeffs(a)
% Object ``a`` must belong to a class derived from
% Transport, and that was constructed using a model that implements
% mixture-averaged transport properties. If not, you will get the error message ::
%
%     **** Method getMixDiffCoeffs not implemented. ****
%
% In this case, create a new gas mixture model that uses a mixture-averaged
% transport manager, for example::
%
%     >> gas = GRI30('mixture-averaged');
%
% See also: :mat:func:`MultiDiffCoeffs`
%
% :param a:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
%     for which mixture-averaged diffusion coefficients are desired.
% :return:
%     Vector of length nSpecies with the mixture-averaged diffusion
%     coefficients. Units: m**2/s
%

v = trans_get(a.id, 11, nSpecies(a.th));
