function v = thermalDiffCoeffs(a)
% THERMALDIFFCOEFFS  Get the thermal diffusion coefficients.
% v = thermalDiffCoeffs(a)
% Object ``a`` must belong to a class derived from
% Transport, and that was constructed by specifying the ``'multicomponent'``
% option. If ``'multicomponent'`` was not specified, the returned values will
% all be zero.
%
% :param a:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
%     for which the thermal diffusion coefficients are desired.
% :return:
%     Vector of thermal diffusion coefficients of length nSpecies
%

v = trans_get(a.id, 12, nSpecies(a.th));
