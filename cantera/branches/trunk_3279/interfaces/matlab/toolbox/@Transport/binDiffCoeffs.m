function v = binDiffCoeffs(a)
% BINDIFFCOEFFS  Get the binary diffusion coefficents.
% v = binDiffCoeffs(a)
% :param a:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
%     for which binary diffusion coefficients are desired.
% :return:
%     A matrix of binary diffusion coefficients.
%     The matrix is symmetric: d(i,j) = d(j,i). Units: m**2/s
%

v = trans_get(a.id, 21, nSpecies(a.th));
