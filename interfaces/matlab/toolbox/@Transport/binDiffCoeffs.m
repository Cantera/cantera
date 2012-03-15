function v = binDiffCoeffs(a)
%BINDIFFCOEFFS  Binary diffusion coefficients (m^2/s).
%
%      d = binDiffCoeffs(gas)
%
%      returns the matrix of binary diffusion coefficients in array
%      d. The matrix is symmetric: d(i,j) = d(j,i).
%
v = trans_get(a.id, 21, nSpecies(a.th));
