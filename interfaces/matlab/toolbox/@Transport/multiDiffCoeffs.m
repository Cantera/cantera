function v = multiDiffCoeffs(a)
% MULTIDIFFCOEFFS  Get the multicomponent diffusion coefficients.
% v = multiDiffCoeffs(a)
% Object ``a`` must belong to a class derived from
% Transport, and that was constructed by specifying the ``'Multi'``
% option. If ``'Multi'`` was not specified, you will get the
% error message ::
%
%     **** Method getMultiDiffCoeffs not implemented. ****
%
% In this case, try method :mat:func:`mixDiffCoeffs`, or create a
% new gas mixture model that uses a mixture-averaged transport manager,
% for example::
%
%     >> gas = GRI30('Multi');
%
% :param a:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
%     for which multicomponent diffusion coefficients are desired.
% :return:
%     Matrix of size [nSpecies, nSpecies] with the multicomponent 
%     diffusion coefficients D_{i,j}. Units: m**2/s
%
%     The coefficient D_{row,column} is used for computing the
%     contribution of concentration gradients in species [column] 
%     to the diffusive flux of species [row].
%

v = trans_get(a.id, 22, nSpecies(a.th));

